# NTC limma results
ann <- read_rds("../Downloads/limma_perCTex.vivovsin.vivo.rds")
ann <- data.frame(ann)
export_by_celltype(
  df = ann,
  output_dir = here("Limma_results_by_celltype_NTC_1"),
  output_file = "Supplementary_Table1.xlsx",
  sheet_columns = c("genes", "logFC", "adj.P.Val", "celltype", "group")
)

#  NTC enrichment
ann <- read_rds("../Downloads/NTC_fgsea.rds")
ann <- data.frame(ann)
ann <- ann %>%
  filter(db == "MSigDB_Hallmark_2020") %>%
  filter(!is.na(NES))
export_by_celltype(
  df = ann,
  output_dir = here("Enrichment_results_NTC_2"),
  output_file = "Supplementary_Table2.xlsx",
  sheet_columns = c("pathway", "padj", "NES", "size", "celltype", "leadingEdge", "db")
)

# KO limma results
ann <- read_rds("../Downloads/limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds")
setDT(ann)
export_by_celltype(
  df = ann,
  output_dir = here("Limma_results_by_celltype_KOs_3"),
  output_file = "Supplementary_Table3.xlsx",
  sheet_columns = c("genotype", "condition", "gene", "logFC", "adj.P.Val", "celltype", "group")
)

# KO enrichment
ann <- fread("../Downloads/fgsea_MSigDB_Hallmark_2020.tsv")
setDT(ann)
ann <- ann[db == "MSigDB_Hallmark_2020"]
ann <- ann[!is.na(NES)]
ann <- data.frame(ann)
export_by_celltype(
  df = ann,
  output_dir = here("Enrichment_results_KOs_4"),
  output_file = "Supplementary_Table4.xlsx",
  sheet_columns = c("pathway", "genotype", "padj", "NES", "size", "celltype", "leadingEdge", "db")
)

# culture effect to KO effect
ann <- read_rds("../Downloads/enrichment_to_NTC_genes.rds")
ann <- data.frame(ann)
ann <- ann %>%
  filter(overlap > 5)
export_by_celltype(
  df = ann,
  output_dir = here("Overlap_in.vivo_KO_to_NTC_5"),
  output_file = "Supplementary_Table5.xlsx",
  sheet_columns = c("coef", "celltype", "overlap", "p.value", "odds.ratio")
)

# Diff_exp_JAK_STAT
ann <- read_rds("../Downloads/combined_jakstat_diff_exp.rds")
ann <- data.frame(ann)
ann <- ann %>%
  mutate(
    genotype = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ str_remove(coef, "^ex\\.vivo_genotype"),
      str_detect(coef, "^genotype.*:treatment") ~ str_remove(str_extract(coef, "genotype[^:]+"), "genotype"),
      str_detect(coef, "^genotype") ~ str_remove(coef, "^genotype"),
      TRUE ~ NA_character_
    ),
    condition = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ "ex.vivo",
      str_detect(coef, ":treatment") ~ "interaction",
      str_detect(coef, "^genotype") ~ "in.vivo",
      TRUE ~ NA_character_
    )
  )
ann$celltype <- ann$cell_type
export_by_celltype(
  df = ann,
  output_dir = here("Diff_exp_JAKSTAT_6"),
  output_file = "Supplementary_Table6.xlsx",
  sheet_columns = c("genotype", "condition", "gene", "logFC", "adj.P.Val", "celltype")
)

# KO enrichment JAKSTAT
ann <- read_rds("../Downloads/fgsea_hom_vs_ex.vivo_per_CT.rds")
setDT(ann)
ann <- ann[db == "MSigDB_Hallmark_2020"]
ann <- ann[!is.na(NES)]
ann <- data.frame(ann)
ann <- ann %>%
  mutate(
    genotype = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ str_remove(coef, "^ex\\.vivo_genotype"),
      str_detect(coef, "^genotype.*:treatment") ~ str_remove(str_extract(coef, "genotype[^:]+"), "genotype"),
      str_detect(coef, "^genotype") ~ str_remove(coef, "^genotype"),
      TRUE ~ NA_character_
    ),
    condition = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ "ex.vivo",
      str_detect(coef, ":treatment") ~ "interaction",
      str_detect(coef, "^genotype") ~ "in.vivo",
      TRUE ~ NA_character_
    )
  )
export_by_celltype(
  df = ann,
  output_dir = here("JAKSTAT_Enrichment_results_7"),
  output_file = "Supplementary_Table7.xlsx",
  sheet_columns = c("pathway", "genotype", "condition", "padj", "NES", "size", "celltype", "leadingEdge", "db")
)
