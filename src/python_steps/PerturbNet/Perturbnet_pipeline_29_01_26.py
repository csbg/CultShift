#!/usr/bin/env python

# ============================================================
# 0. Imports
# ============================================================

import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

import scanpy as sc
import anndata as ad
import scvi

import torch
import torch.nn as nn
import torch.nn.functional as F

import mygene
import matplotlib.pyplot as plt

from perturbnet.util import *
from perturbnet.cinn.flow import *
from perturbnet.genotypevae.genotypeVAE import *
from perturbnet.data_vae.vae import *
from perturbnet.cinn.flow_generate import SCVIZ_CheckNet2Net


# ============================================================
# 1. Load AnnData (RAW COUNTS)
# ============================================================

adata = ad.read_h5ad("data/TfCf_All_rawcounts.h5ad")

print("Loaded AnnData:", adata)
print("X dtype:", adata.X.dtype)
print("X max:", adata.X.max())


# ============================================================
# 2. Load and align metadata
# ============================================================

import pandas as pd

# Load metadata
meta_ex = pd.read_csv("data/ex.vivo_metadata.tsv", sep="\t")
meta_in = pd.read_csv("data/in.vivo_metadata.tsv", sep="\t")

print(meta_ex.head())
print(meta_ex.columns)

# Concatenate ex vivo and in vivo metadata
meta = pd.concat([meta_ex, meta_in], axis=0, ignore_index=True)

assert meta["rn"].is_unique

# Extract core 10x barcode
adata.obs["barcode"] = adata.obs_names.str[-18:]  # last 18 chars
meta["barcode"] = meta["rn"].str[:18]            # first 18 chars

# Check how many are common
common = adata.obs["barcode"].isin(meta["barcode"])
print(f"Common cells: {common.sum()} / {adata.n_obs}")

# Keep only metadata rows that exist in adata
meta_sub = meta[meta["barcode"].isin(adata.obs["barcode"])].copy()

# Remove duplicates in metadata
meta_sub = meta_sub.drop_duplicates(subset="barcode")

# Set barcode as index for alignment
meta_sub = meta_sub.set_index("barcode")

# Subset adata.obs to only overlapping barcodes
adata = adata[adata.obs["barcode"].isin(meta_sub.index)].copy()

# Reorder metadata to match adata
meta_sub = meta_sub.loc[adata.obs["barcode"]]

print("AnnData shape:", adata.shape)
print("Metadata shape:", meta_sub.shape)

# Confirm all barcodes aligned
all_match = (adata.obs["barcode"] == meta_sub.index).all()
print("All barcodes aligned:", all_match)

# ============================================================
# 3. Attach metadata to AnnData
# ============================================================
# ============================================================
# 3. Attach metadata to AnnData
# ============================================================

adata.obs["guide"] = meta_sub["guide"].values
adata.obs["cell_type"] = meta_sub["functional.cluster"].values
adata.obs["tissue"] = meta_sub["tissue"].values
adata.obs["timepoint"] = meta_sub["timepoint"].values
adata.obs["mixscape_class"] = meta_sub["mixscape_class.global"].values

adata.obs["genotype"] = (
    adata.obs["guide"].str.split("_").str[0].astype("category")
)
# 4. GO onehot matrix
#============================================================
adata.obs["genotype"] = adata.obs["genotype"].astype("category")

perturb_genes = sorted(
    g for g in adata.obs['genotype'].cat.categories
    if g not in ['NTC', 'NA']  # remove controls or missing
)
import mygene

# ======================================
# 3. Query GO annotations using mygene
# ======================================
mg = mygene.MyGeneInfo()
res = mg.querymany(
    perturb_genes,
    scopes='symbol',
    fields='go',
    species='mouse'  # change if human
)

# Build gene → GO dictionary
gene_to_go = {}
for r in res:
    gene = r.get('query')
    if 'go' not in r:
        gene_to_go[gene] = []
        continue
    gos = []
    for domain in ['BP', 'MF', 'CC']:
        if domain in r['go']:
            entries = r['go'][domain]
            if isinstance(entries, list):
                gos.extend([e['id'] for e in entries])
            else:
                gos.append(entries['id'])
    gene_to_go[gene] = list(set(gos))  # remove duplicates

# Remove genes with zero GO annotations
perturb_genes = [g for g in perturb_genes if sum([1 for go in gene_to_go[g]]) > 0]
# Rebuild GO vocabulary
all_go_terms = sorted({go for gos in gene_to_go.values() for go in gos})

gene_to_idx = {gene: i for i, gene in enumerate(perturb_genes)}
go_to_idx   = {go: i for i, go in enumerate(all_go_terms)}

onehot_matrix = np.zeros(
    (len(perturb_genes), len(all_go_terms)), dtype=np.int8
)

for gene, gos in gene_to_go.items():
    if gene not in gene_to_idx:
        continue
    g_idx = gene_to_idx[gene]
    for go in gos:
        onehot_matrix[g_idx, go_to_idx[go]] = 1

# mostly zeros, ~10–200 ones
i = 0
print("Gene:", perturb_genes[i])
print("Number of GO terms:", onehot_matrix[i].sum())
go_indices = np.where(onehot_matrix[i] == 1)[0]
go_terms = [all_go_terms[j] for j in go_indices[:10]]
print(go_terms)
zero_go_genes = [
    perturb_genes[i]
    for i in range(len(perturb_genes))
    if onehot_matrix[i].sum() == 0
]
print(zero_go_genes)
#testplot
import matplotlib.pyplot as plt
import os

# Output directory
outdir = "results/plots"
os.makedirs(outdir, exist_ok=True)

plt.figure(figsize=(6, 4))
plt.hist(onehot_matrix.sum(axis=1), bins=30)
plt.xlabel("GO terms per gene")
plt.ylabel("Number of genes")
plt.title("Distribution of GO terms per perturbation gene")

plt.tight_layout()
plt.savefig(os.path.join(outdir, "go_terms_per_gene_hist.png"), dpi=300)
plt.show()
#======================================================
# 5. Ensure raw counts (INTEGER) + counts layer
# ============================================================

adata.X = sp.csr_matrix(adata.X.astype(np.int32))
adata.layers["counts"] = adata.X.copy()

print("Counts layer dtype:", adata.layers["counts"].dtype)


# ============================================================
# 6. Train / test split
# ============================================================

train_mask = (
    ((adata.obs["tissue"] == "ex.vivo") & (adata.obs["genotype"] != "NA")) |
    ((adata.obs["tissue"] == "ex.vivo") & (adata.obs["genotype"] == "NTC")) |
    ((adata.obs["tissue"] == "in.vivo") & (adata.obs["genotype"] == "NTC"))
)

test_mask = (
    (adata.obs["tissue"] == "in.vivo") &
    (adata.obs["genotype"] != "NTC") &
    (adata.obs["genotype"] != "NA")
)

adata_train = adata[train_mask].copy()
adata_test  = adata[test_mask].copy()

print("Train cells:", adata_train.n_obs)
print("Test cells :", adata_test.n_obs)


# ============================================================
# 7. scVI setup and training
# ============================================================

scvi.data.setup_anndata(
    adata_train,
    layer="counts"
)

scvi_model = scvi.model.SCVI(
    adata_train,
    n_latent=10
)

scvi_model.train(
    n_epochs=500,
    frequency=20
)

os.makedirs("results/model", exist_ok=True)
scvi_model.save("results/model")

print("✔ scVI model trained and saved")


# ============================================================
# 8. Plot training loss
# ============================================================

train_loss = np.array(scvi_model.history["elbo_train_set"])
valid_loss = np.array(scvi_model.history["elbo_test_set"])

epochs = np.arange(len(train_loss)) * 20

plt.figure(figsize=(8, 6))
plt.plot(epochs, train_loss, label="Train")
plt.plot(epochs, valid_loss, label="Validation")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("scVI Training Loss")
plt.legend()
plt.tight_layout()
plt.show()


# ============================================================
# DONE
# ============================================================

print("✅ Pipeline completed successfully")
