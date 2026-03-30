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

basedir <- dirout("Ag_ScRNA_15_celltype_biolord_limma_v2_mye/")

#####################################################################
#load data and clean meta_predicted
################################################
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
#meta_predicted
meta <- fread(inDir1("meta_cleaned.tsv"))
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE]
grep("Brd9_BR_48005 MkP inVivo_OP1_lin-_14d_1", rownames(meta), value = T)
#ex.vivo from opriginal meta
meta_ex <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "ex.vivo")
#in vivo NTC from original meta
meta_in_NTC <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "in.vivo", genotype == "NTC")
#meta with in vivo and ex vivo NTC + ex vivo KO
meta <- rbind(meta_ex,meta_in_NTC)

in.vivo_predicted <- as.data.frame(counts(read_rds("/vscratch/aarathy/tfcf-biolord/data_generated/model_guide/predicted_sum.rds")))
# Find duplicated column names
# Ensure in.vivo_predicted is a data.frame
in.vivo_predicted <- as.data.frame(in.vivo_predicted)

# Collapse duplicates by taking the mean
in.vivo_predicted <- sapply(
  unique(colnames(in.vivo_predicted)),
  function(x) {
    tmp <- in.vivo_predicted[, colnames(in.vivo_predicted) == x, drop = FALSE]  # always a matrix
    rowMeans(tmp)
  }
)

# Convert back to matrix
in.vivo_predicted <- as.matrix(in.vivo_predicted)

# Check
any(duplicated(colnames(in.vivo_predicted)))  # should be FALSE
dim(in.vivo_predicted)

# Check that all column names are now unique
any(duplicated(colnames(in.vivo_predicted)))

meta_predicted <- read_rds("/vscratch/aarathy/tfcf-biolord/data_generated/model_guide/predicted_sum.rds")
meta_predicted <- as.data.frame(colData(meta_predicted))
meta_predicted <- meta_predicted %>%
  mutate(celltype = celltype_projection)


# Modify celltype column in meta_predicted
meta_predicted$celltype <- gsub("Gran. P", "Gran.P", meta_predicted$celltype)
meta_predicted$celltype <- gsub("Eo/Ba", "Eo.Ba", meta_predicted$celltype)

# Modify row names of meta_predicted
rownames(meta_predicted) <- gsub("Gran. P", "Gran.P", rownames(meta_predicted))
rownames(meta_predicted) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta_predicted))


# Modify column names of in.vivo_predicted
colnames(in.vivo_predicted) <- gsub("Gran. P", "Gran.P", colnames(in.vivo_predicted))
colnames(in.vivo_predicted) <- gsub("Eo/Ba", "Eo.Ba", colnames(in.vivo_predicted))

#colnames(in.vivo_predicted) <- sub("\\.[0-9]+$", "", colnames(in.vivo_predicted))
in.vivo_predicted <- in.vivo_predicted[,unique(colnames(in.vivo_predicted))]
colnames(in.vivo_predicted) <- make.unique(colnames(in.vivo_predicted))

all(colnames(in.vivo_predicted) %in% row.names(meta_predicted))
all(row.names(meta_predicted) %in% colnames(in.vivo_predicted))
setdiff(row.names(meta_predicted), colnames(in.vivo_predicted))


library(Matrix)
library(dplyr)

# collapse duplicated columns by summing counts
in.vivo_predicted <- in.vivo_predicted[, !duplicated(colnames(in.vivo_predicted))]  # simple remove duplicates




meta_predicted <- meta_predicted[rownames(meta_predicted)[rownames(meta_predicted) %in% colnames(in.vivo_predicted)],]
#exclude NTC from the predicted set metadata and data
meta_predicted <- meta_predicted %>%
  filter(genotype != "NTC")


meta_predicted <- meta_predicted %>%
  filter(celltype != "Imm. B-cell") %>%
  dplyr::select(-c(celltype_projection))

meta_predicted <- meta_predicted %>%
  dplyr::rename("tissue" = "dataset")%>%
  mutate(genotype = gsub("_.*$","",guide))%>%
  dplyr::select(genotype,tissue,celltype)

in.vivo_predicted <- in.vivo_predicted[, rownames(meta_predicted)]
#bind observed(exvivoNTC, in vivo NTC, ex vivo KO) and predicted (in vivo KO) metadata
meta <- rbind(meta,meta_predicted) %>%
  filter(genotype %in% c("NTC", unique(meta_predicted$genotype))) %>%
  filter(celltype %in% unique(meta_predicted$celltype))

meta$tissue <- gsub("invivo", "in.vivo", meta$tissue)
#counts
rownames(meta[meta$tissue=="ex.vivo",])
rownames(meta[meta$tissue=="in.vivo",])
colnames(in.vivo_predicted)
counts <- read.delim(inDir("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
grep("Brd9_BR_48005.MkP.inVivo_OP1_lin._14d_1", colnames(counts), value = T)
Smarcd1_AS_15166_inVivo_OP4_lin-_14d_1
grep("Brd9_BR_48005 MkP inVivo_OP1_lin-_14d_1", rownames(meta), value = T)
exvivo <- rownames(meta[meta$tissue == "ex.vivo",])
invivo_NTC <-  meta %>%
  filter(genotype =="NTC")%>%
  filter(tissue == "in.vivo")%>%
  rownames()
subset <- c(exvivo, invivo_NTC)
counts <- counts[!(rownames(counts) %in% genes_to_exclude),subset]

counts <- cbind(counts,in.vivo_predicted[rownames(counts),rownames(meta_predicted)]) 
stopifnot(all(colnames(counts)==rownames(meta)))

stopifnot(setequal(unique(meta[meta$tissue == "ex.vivo", ]$celltype), 
                   unique(meta[meta$tissue == "in.vivo", ]$celltype)))

################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
meta%>% write_rds(basedir("meta.rds"))
################################################
#perform DE independent of celltype
################################################
#out <- dirout(paste0(base,"/per_cell_type/"))
model_formula <- "~ tissue * genotype"
limmaRes <- performDE(meta, counts, model_formula)

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

limmaRes %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred.rds"))
limmaRes_int <- limmaRes[limmaRes$coef %in% grep("interaction",limmaRes$coef,value=T),]%>%na.omit()
#dataVoom<-result$dataVoom
limmaRes_int %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_interaction_pred.rds"))
###############################################################################  
unique(limmaRes$celltype)
