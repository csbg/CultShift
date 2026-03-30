####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#limma with all 
#data
#
#Aarathy
###############
###############
source("src/00_init.R")
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(stringr)
#####################################################################
inDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
inDir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
inDir2 <- dirout("Ag_ScRNA_23_PerturbNet_python_steps")
basedir <- dirout("Ag_ScRNA_26_prediction_COMBINED_dataVoom/")

#####################################################################
#load data and clean meta_biolord
################################################
#in.vivo data---------
counts <- read.delim(inDir("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
in.vivo_biolord <- as.data.frame(counts(read_rds("/vscratch/aarathy/tfcf-biolord/data_generated/model_guide/predicted_sum.rds")))
pred <- read.csv(
  inDir2("predicted_expression_SUM_direct_scaled_to_invivo_epoch35.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = NULL
)
str(pred)

pred$celltype <- pred$cell_type
pred$cell_type <- NULL
# Step 2: Combine guide, sample, cell_type, tissue to make rownames
pred$combined_name <- paste0(pred$guide,".", pred$celltype,".", pred$sample)
# Find duplicated combined names
dup_names <- pred$combined_name[duplicated(pred$combined_name)]
# Step 3: Assign combined_name as rownames
rownames(pred) <- pred$combined_name
non_numeric_cols <- !sapply(pred, is.numeric)

# Names of non-numeric columns
names(pred)[non_numeric_cols]
# Find the column index of "Mrpl15"
# counts_table
perturb_counts <- pred[, !non_numeric_cols]
perturb_counts <- t(perturb_counts)


##############
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
#meta_invivo
meta_obs <- fread(inDir1("meta_cleaned.tsv"))
meta_obs <- as.data.frame(meta_obs)               # Convert to dataframe (optional)
rownames(meta_obs) <- meta_obs[[1]]   

meta_obs <- meta_obs[, -1, drop = FALSE]

#in vivo NTC from original meta
meta <- meta_obs %>%
  dplyr::select(genotype,tissue,celltype) 
meta$model <- "Observed"

counts <- counts[!(rownames(counts) %in% genes_to_exclude),rownames(meta)]
#biolord------------------------

meta_biolord <- read_rds("/vscratch/aarathy/tfcf-biolord/data_generated/model_guide/predicted_sum.rds")
meta_biolord <- as.data.frame(colData(meta_biolord))
meta_biolord <- meta_biolord %>%
  mutate(celltype = celltype_projection)

# ensures the indices match and then make unique names for in.vivo predicted(
#actuaklly they were different samples. since the order is preserved proceeded like below to make colnames unique

all(
  make.unique(colnames(in.vivo_biolord)) ==
    make.unique(rownames(meta_biolord))
)

colnames(in.vivo_biolord) <- make.unique(colnames(in.vivo_biolord))


# Modify celltype column in meta_biolord
meta_biolord$celltype <- gsub("Gran. P", "Gran.P", meta_biolord$celltype)
meta_biolord$celltype <- gsub("Eo/Ba", "Eo.Ba", meta_biolord$celltype)

# Modify row names of meta_biolord
rownames(meta_biolord) <- gsub("Gran. P", "Gran.P", rownames(meta_biolord))
rownames(meta_biolord) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta_biolord))


# Modify column names of in.vivo_biolord
colnames(in.vivo_biolord) <- gsub("Gran. P", "Gran.P", colnames(in.vivo_biolord))
colnames(in.vivo_biolord) <- gsub("Eo/Ba", "Eo.Ba", colnames(in.vivo_biolord))


in.vivo_biolord <- in.vivo_biolord[,unique(colnames(in.vivo_biolord))]

all(colnames(in.vivo_biolord) %in% row.names(meta_biolord))
all(row.names(meta_biolord) %in% colnames(in.vivo_biolord))

meta_biolord <- meta_biolord[rownames(meta_biolord)[rownames(meta_biolord) %in% colnames(in.vivo_biolord)],]
#exclude NTC from the predicted set metadata and data
meta_biolord <- meta_biolord %>%
  filter(genotype != "NTC") 

meta_biolord <- meta_biolord %>%
  filter(celltype != "Imm. B-cell") %>%
  dplyr::select(-c(celltype_projection))

meta_biolord <- meta_biolord %>%
  dplyr::rename("tissue" = "dataset")%>%
  mutate(genotype = gsub("_.*$","",guide))%>%
  dplyr::select(genotype,tissue,celltype)
meta_biolord$model <- "biolord"
############



in.vivo_biolord <- in.vivo_biolord[, rownames(meta_biolord)]

stopifnot(identical(rownames(meta_biolord), colnames(in.vivo_biolord)))

#bind observed(exvivoNTC, in vivo NTC, ex vivo KO) and predicted (in vivo KO) metadata

##############################################################################
#PerturbNet-------------------------------------------
# Correct way to read your CSV (all columns are regular columns)


# Step 4: Remove the old metadata columns (including combined_name)
#build metadata for pred------
meta_perturb <- data.frame(
  sample = rownames(pred)
)

meta_perturb$genotype <- pred$genotype
meta_perturb$celltype <- pred$celltype
meta_perturb$tissue   <- pred$tissue
rownames(meta_perturb) <- meta_perturb$sample
meta_perturb$rowname <- rownames(meta_perturb)

# Step 5: Confirm column names and structure
# Logical vector: TRUE if column is NOT numeric


meta_perturb <- meta_perturb %>%
  mutate(celltype = celltype %>%
           gsub("^MEP \\(G1\\)$", "MEP.G1", .) %>%
           gsub("^MEP \\(S\\)$", "MEP.S", .) %>%
           gsub("^MEP \\(early\\)$", "MEP.early", .) %>%
           gsub("^MEP \\(pert\\.\\)$", "MEP.pert.", .) %>%
           gsub("^GMP \\(early\\)$", "GMP.early", .) %>%
           gsub("^GMP \\(late\\)$", "GMP.late", .) %>%
           gsub("^Gran\\. P$", "Gran.P", .) %>%
           gsub("^Imm\\. B-cell$", "Imm.B.cell", .) %>%
           gsub("^Eo\\/Ba$", "Eo.Ba", .)
         
  )



all(rownames(meta_perturb) == colnames(perturb_counts))
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell",
                          "Imm.B.cell","MEP","MEP.pert.","MEP.S", "MEP.G1","GMP.late",
                          "GMP.early","Imm. B-cell","nan")

meta_perturb <- meta_perturb%>%
  filter(!(celltype %in% celltypes_to_exclude))
meta_perturb <- meta_perturb %>%
  filter(tissue == "in.vivo") %>%
  filter(!(genotype == "NTC"))

perturb_counts <- perturb_counts[,rownames(meta_perturb)]
meta_perturb$rowname <- NULL
meta_perturb$sample <- NULL
#important to match real in vivo meta data
rownames(meta_perturb) <- gsub("[\\ \\(\\)-]", ".", rownames(meta_perturb))
colnames(perturb_counts) <-  gsub("[\\ \\(\\)-]", ".", colnames(perturb_counts))
rownames(meta_perturb) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta_perturb))
colnames(perturb_counts) <-  gsub("Eo/Ba", "Eo.Ba", colnames(perturb_counts))


#observed metadata

#ex.vivo from opriginal meta
meta_ex <- meta_obs %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "ex.vivo")

meta_in_NTC <- meta_obs %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "in.vivo", genotype == "NTC")

meta_in_KO <- rownames(meta_obs)[!(rownames(meta_obs) %in% c(rownames(meta_in_NTC), rownames(meta_ex)))]
table(rownames(meta_perturb) %in% meta_in_KO)

meta_perturb <- meta_perturb[rownames(meta_perturb) %in% meta_in_KO,]
meta_perturb$model <-"PerturbNet"
perturb_counts <- perturb_counts[,rownames(meta_perturb)]
#data
meta_perturb <- meta_perturb%>%
  filter(genotype %in% unique(meta_ex$genotype))
meta_perturb <- meta_perturb%>%
  filter(celltype %in% unique(meta_ex$celltype))

# PerturbNet
colnames(perturb_counts) <- paste0("PerturbNet_", colnames(perturb_counts))
rownames(meta_perturb)   <- paste0("PerturbNet_", rownames(meta_perturb))

#combine
stopifnot(all(rownames(meta_biolord) == colnames(in.vivo_biolord)))
stopifnot(all(rownames(meta) == colnames(counts)))
stopifnot(all(rownames(meta_perturb) == colnames(perturb_counts)))

meta_all <- rbind(meta,meta_biolord) 
meta_all$tissue <- gsub("invivo", "in.vivo", meta$tissue)
meta_all <- rbind(meta_all,meta_perturb)

meta_all %>% write_rds(basedir("meta_all.rds"))


# subsetiing counts
# Ensure same genes and order
common_genes <- Reduce(
  intersect,
  list(
    rownames(counts),           # Observed
    rownames(in.vivo_biolord),  # biolord
    rownames(perturb_counts)    # PerturbNet
  )
)

counts_all <- cbind(
  counts[common_genes, ],
  in.vivo_biolord[common_genes, ],
  perturb_counts[common_genes, ]
)


counts_all <- counts_all[, rownames(meta_all)]

stopifnot(all(rownames(meta_all) == colnames(counts_all)))
#
meta_all$group <- with(
  meta_all,
  paste(celltype, model, genotype, tissue, sep = "_")
)

meta_all$group <- factor(meta_all$group)

design <- model.matrix(~ 0 + group, data = meta_all)
colnames(design) <- levels(meta_all$group)

dge <- DGEList(counts_all)
dge <- calcNormFactors(dge)

dataVoom <- voom(dge, design, plot = FALSE)
stopifnot(all(rownames(meta_all) == colnames(counts_all)))
dataVoom %>% write_rds(basedir("dataVoom.rds"))
meta_all %>% write_rds(basedir("meta_all.rds"))

#######################


contrast_list <- list()

for (ct in unique(meta_all$celltype)) {
  
  ref <- paste(ct, "Observed", "NTC", sep = "_")
  
  if (!ref %in% colnames(design)) next
  
  test_groups <- grep(
    paste0("^", ct, "_"),
    colnames(design),
    value = TRUE
  )
  
  test_groups <- setdiff(test_groups, ref)
  
  for (tg in test_groups) {
    cname <- paste(tg, "vs", ref, sep = "_")
    contrast_list[[cname]] <- paste(tg, "-", ref)
  }
}

contrast_matrix <- makeContrasts(
  contrasts = unlist(contrast_list),
  levels = design
)

dge <- DGEList(counts_all)
dge <- calcNormFactors(dge)

dataVoom <- voom(dge, design, plot = FALSE)
dataVoom %>% write_rds(basedir("dataVoom.rds"))
