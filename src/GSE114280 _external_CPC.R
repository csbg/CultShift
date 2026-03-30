source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(tibble)
library(data.table)
library(msigdbr)
library(org.Hs.eg.db) 
library(fgsea)# for human; use org.Mm.eg.db for mouse

InDir1 <- dirout("Ag_top_filtered_genes")

out <- dirout("GSE114280_CPC_external_dataset")
fresh_CDC <- fread(out("../../GSE114280_CPC/GSM3138748_fresh_CPC.txt/GSM3138748_fresh_CPC.txt"))
cultured_CDC <- fread(out("../../GSE114280_CPC/GSM3138749_cultured_CPC.txt/GSM3138749_cultured_CPC.txt"))
str(fresh_CDC)
# Create pseudobulk vectors
pseudobulk1 <- rowSums(fresh_CDC[, -1, with = FALSE])
names(pseudobulk1) <- fresh_CDC$V1

pseudobulk2 <- rowSums(cultured_CDC[, -1, with = FALSE])
names(pseudobulk2) <- cultured_CDC$V1

# Combine into a dataframe
# First, make sure gene order matches
all_genes <- union(names(pseudobulk1), names(pseudobulk2))

df_pseudobulk <- data.frame(
  gene = all_genes,
  fresh_CDC = pseudobulk1[all_genes],
  cultured_CDC = pseudobulk2[all_genes]
)

# View
head(df_pseudobulk)
rownames(df_pseudobulk) <- df_pseudobulk$gene

#change ID to names

# convert counts to data.table if not already
counts <- data.frame(df_pseudobulk)%>%na.omit()

metadata <- data.frame(
  sample = colnames(counts)
) %>%
  mutate(
    tissue = ifelse(grepl("cul", sample), "ex.vivo",
                    "in.vivo"))
  

rownames(metadata) <- metadata$sample

# Make sure factor levels match actual values
metadata$tissue <- factor(metadata$tissue, 
                          levels = c("in.vivo","ex.vivo"))

# Double check order
stopifnot(all(colnames(counts) == rownames(metadata)))

# Create design matrix
design <- model.matrix(~tissue, data = metadata)
rownames(design) <- rownames(metadata)
rownames(counts) <- counts$gene
counts <- counts[,c("fresh_CDC", "cultured_CDC")]
counts <- as.matrix(counts) 

##################################
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d <- d0[keep,]


###############################################################################
#setting the model
###############################################################################


# Normalization
dataVoom <- voom(d, design, plot = T)
dataVoom %>% write_rds(out("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

#############
genes_fig1 <- read_rds(InDir1("genes_fig1.rds"))
head(genes_fig1)
unique(genes_fig1$pathway)
ifn_genes <- genes_fig1 %>%
  filter(pathway =="ISG core")%>%
  pull(genes)

sign_genes <- intersect(ifn_genes, rownames(dataVoom$E))


dat.list <- list()
for(gg in sign_genes){
  dat.list[[gg]] <- metadata %>%
    mutate(E=dataVoom$E[gg,]) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}


sign <- bind_rows(dat.list, .id ="genes")
  #left_join(sign_genes_tbl, by = "genes")  # keep coef info


mat <- sign %>%
  dplyr::select(sample, genes, E) %>%
  pivot_wider(
    names_from = sample,
    values_from = E,
    values_fn = mean  # <- take mean if duplicates exist
  ) %>%
  column_to_rownames("genes") %>%
  as.matrix()



hc_genes <- hclust(dist(mat))
gene_order <- hc_genes$labels[hc_genes$order]
sign$genes <- factor(sign$genes, levels = gene_order)

sign$tissue <- factor(sign$tissue, levels = c(
  "in.vivo","ex.vivo"
))
sign$sample <- factor(sign$sample, levels = c(
  "fresh_CDC","cultured_CDC"
))
(sign_plot <- ggplot(sign, aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
   # facet_grid(~ tissue, space="free", scales="free") +
    scale_fill_gradient2(low="blue", high="red") +
    optimized_theme_fig())

###########
#mtorc1
genes_fig1 <- read_rds(InDir1("genes_fig1.rds"))
head(genes_fig1)
unique(genes_fig1$pathway)
chol_genes <- genes_fig1 %>%
  filter(pathway %in% c("mTORC1_or_Cholesterol",
                        "mTORC1/Cholesterol" ))%>%
  pull(genes)

sign_genes <- intersect(toupper(chol_genes), rownames(dataVoom$E))
dat.list <- list()
for(gg in sign_genes){
  dat.list[[gg]] <- metadata %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
sign_genes_tbl <- limmaRes %>%
  filter(genes %in% sign_genes)

sign <- bind_rows(dat.list, .id ="genes") %>%
  left_join(sign_genes_tbl, by = "genes")  # keep coef info


mat <- sign %>%
  dplyr::select(sample, genes, E) %>%
  pivot_wider(
    names_from = sample,
    values_from = E,
    values_fn = mean  # <- take mean if duplicates exist
  ) %>%
  column_to_rownames("genes") %>%
  as.matrix()



hc_genes <- hclust(dist(mat))
gene_order <- hc_genes$labels[hc_genes$order]
sign$genes <- factor(sign$genes, levels = gene_order)


(sign_plot <- ggplot(sign, aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(coef ~ tissue, space="free", scales="free") +
    scale_fill_gradient2(low="blue", high="red") +
    optimized_theme_fig())
