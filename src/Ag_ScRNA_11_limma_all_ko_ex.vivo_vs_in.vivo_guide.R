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
#####################################################################
InDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
InDir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")

basedir <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_guide/")

#####################################################################
#load data and clean metadata
################################################
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
#metadata
meta <- fread(InDir1("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   
meta <- meta[, -1, drop = FALSE] 
#colnames(meta) <- gsub("rowname","sample1", colnames(meta))
meta$sample1 <- rownames(meta)
# Check if there are at least 2 distinct samples per tissue for each genotype and celltype
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)

replicates_per_ko <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(
    valid_ko = any(valid_ko),
    total_in_vivo = sum(in.vivo, na.rm = TRUE),
    total_ex_vivo = sum(ex.vivo, na.rm = TRUE)
    , .groups = "drop") %>%
  mutate(coef = genotype)# replicates_per_ko <- meta %>%

meta <- meta %>%
  inner_join(ko_flags, by =c ("genotype","celltype"))%>%
  filter(valid_ko)
rownames(meta) <- meta$sample1
#counts
counts <- read.delim(InDir("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
counts <- counts[!(rownames(counts) %in% genes_to_exclude),rownames(meta)]

stopifnot(all(colnames(counts)==rownames(meta)))

stopifnot(setequal(unique(meta[meta$tissue == "ex.vivo", ]$celltype), 
                   unique(meta[meta$tissue == "in.vivo", ]$celltype)))

################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
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

limmaRes %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
limmaRes_int <- limmaRes[limmaRes$coef %in% grep("interaction",limmaRes$coef,value=T),]%>%na.omit()
#dataVoom<-result$dataVoom
limmaRes_int %>% write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
###############################################################################  

