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
library(data.table)
basedir <- dirout("Ag_ScRNA_23_PerturbNet")
inDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
inDir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
inDir2 <- dirout("Ag_ScRNA_23_PerturbNet_python_steps")


# Correct way to read your CSV (all columns are regular columns)
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


# Step 4: Remove the old metadata columns (including combined_name)
#build metadata for pred------
meta_pred <- data.frame(
  sample = rownames(pred)
)

meta_pred$genotype <- pred$genotype
meta_pred$celltype <- pred$celltype
meta_pred$tissue   <- pred$tissue
rownames(meta_pred) <- meta_pred$sample
meta_pred$rowname <- rownames(meta_pred)


# Step 5: Confirm column names and structure
# Logical vector: TRUE if column is NOT numeric
non_numeric_cols <- !sapply(pred, is.numeric)

# Names of non-numeric columns
names(pred)[non_numeric_cols]
# Find the column index of "Mrpl15"


# counts_table
pred_counts <- pred[, !non_numeric_cols]
pred_counts <- t(pred_counts)




meta_pred <- meta_pred %>%
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



all(rownames(meta_pred) == colnames(pred_counts))
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell",
                          "Imm.B.cell","MEP","MEP.pert.","MEP.S", "MEP.G1","GMP.late",
                          "GMP.early","Imm. B-cell","nan")

meta_pred <- meta_pred%>%
  filter(!(celltype %in% celltypes_to_exclude))
meta_pred <- meta_pred %>%
  filter(tissue == "in.vivo") %>%
  filter(!(genotype == "NTC"))

pred_counts <- pred_counts[,rownames(meta_pred)]
meta_pred$rowname <- NULL
meta_pred$sample <- NULL
#important to match real in vivo meta data
rownames(meta_pred) <- gsub("[\\ \\(\\)-]", ".", rownames(meta_pred))
colnames(pred_counts) <-  gsub("[\\ \\(\\)-]", ".", colnames(pred_counts))
rownames(meta_pred) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta_pred))
colnames(pred_counts) <-  gsub("Eo/Ba", "Eo.Ba", colnames(pred_counts))


#observed metadata
meta_obs  <- fread(inDir1("meta_cleaned.tsv"))
meta_obs <- as.data.frame(meta_obs)               # Convert to dataframe (optional)
rownames(meta_obs) <- meta_obs[[1]]   
meta <- meta_obs
meta <- meta[, -1, drop = FALSE]
#ex.vivo from opriginal meta
meta_ex <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "ex.vivo")

meta_in_NTC <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "in.vivo", genotype == "NTC")

meta <- rbind(meta_ex,meta_in_NTC)

meta_in_KO <- rownames(meta_obs)[!(rownames(meta_obs) %in% rownames(meta))]
table(rownames(meta_pred) %in% meta_in_KO)

meta_pred <- meta_pred[rownames(meta_pred) %in% meta_in_KO,]
pred_counts <- pred_counts[,rownames(meta_pred)]
#data
counts <- read.delim(inDir("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
counts <- counts[,rownames(meta)]


t <- meta[meta$tissue == "ex.vivo",]


meta_pred <- meta_pred%>%
  filter(genotype %in% unique(t$genotype))
meta_pred <- meta_pred%>%
  filter(celltype %in% unique(t$celltype))
#combine
meta_pred <- meta_pred[,colnames(meta)]
meta_all <- rbind(meta,meta_pred)

meta_all <- meta_all %>%
  filter(celltype %in% unique(t$celltype))
#filter data

common_genes <- intersect(rownames(counts), rownames(pred_counts))

counts_all <- cbind(
  counts[common_genes, rownames(meta)],
  pred_counts[common_genes, rownames(meta_pred)]
)

counts_all <- counts_all[,rownames(meta_all)]
all(rownames(meta_all) == colnames(counts_all))
unique(meta_all$celltype)

###############
#baselines
meta_all$genotype <- factor(meta_all$genotype, levels=c("NTC", 
                                                        unique(setdiff(meta_all$genotype,"NTC"))))
meta_all$tissue <- factor(meta_all$tissue, levels=c("in.vivo", "ex.vivo"))
meta_all%>% write_rds(basedir("meta_obs_pred.rds"))


###############
#limma-------------

design <- model.matrix(~ tissue * genotype, data = meta_all)
model_formula <- "~ tissue * genotype"
limmaRes <- performDE(meta = meta_all, counts = counts_all, model_formula = model_formula, thres = 1)
limmaRes <-limmaRes%>%
  filter(coef != "(Intercept)")

limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "^genotype", "in.vivo")) %>%
  mutate(coef = str_replace(coef, "tissue", ""))

limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "ex.vivo.genotype", "interaction"))
# Modify the coefficients
limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    # Replace "interactionX + genotypeX" with "ex.vivoX"
    str_detect(coef, "interaction") & str_detect(coef, "\\+ genotype") ~ 
      str_replace(coef, "interaction", "ex.vivo") %>% 
      str_replace(" \\+ genotype.*", ""),  # Escape + to treat it as a literal character
    # Keep other coefficients unchanged
    TRUE ~ coef
  ))

limmaRes %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred_perturb.rds"))
limmaRes_int <- limmaRes[limmaRes$coef %in% grep("interaction",limmaRes$coef,value=T),]%>%na.omit()
limmaRes_int %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_interaction_pred_perturb.rds"))
unique(limmaRes$celltype)

