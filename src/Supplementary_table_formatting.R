library(dplyr)
library(data.table)
library(openxlsx)
library(here)
library(readr)
library(stringr)

# --- Formatting functions ---
format_padj <- function(padj){
  ret <- ifelse(is.na(padj), NA, format(padj, scientific = TRUE, digits = 3))
  ret <- gsub("e", "E", ret)
  ret <- gsub("^ ", "", ret)
  ret
}

format_numbers <- function(df){
  for(cx in seq_len(ncol(df))){
    if(is.numeric(df[[cx]])){
      df[[cx]] <- round(df[[cx]], 3)
    }
  }
  df
}
# FUNCTIONS ---------------------------------------------------------------
dirout <- function(out, ext="", init=TRUE){
  out.dir <- paste0(here(), "/", out, "/")  # here() not here
  if(init){
    message("Using output directory: ", out.dir)  # no dir.create
  }
  function(...){
    paste0(out.dir, paste0(...), ext)
  }
}



# --- General Excel export function ---
export_by_celltype <- function(df, 
                               output_dir, 
                               output_file, 
                               sheet_columns = NULL, 
                               freeze_first_row = TRUE) {
  
  # Select only relevant columns if provided
  if(!is.null(sheet_columns)){
    df <- df[, intersect(sheet_columns, colnames(df)), drop = FALSE]
  }
  
  df <- as.data.frame(df)
  cell_types <- unique(df$celltype)
  
  wb <- createWorkbook()
  
  for(ct in cell_types){
    ann_ct <- df %>% filter(celltype == ct)
    
    # Format numeric columns
    if("adj.P.Val" %in% colnames(ann_ct)){
      ann_ct$adj.P.Val <- format_padj(ann_ct$adj.P.Val)
    }
    ann_ct <- format_numbers(ann_ct)
    
    # Safe sheet name
    sheet_name <- substr(gsub("[\\/:*?\\[\\]]", "_", ct), 1, 31)
    
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, ann_ct, rowNames = FALSE)
    
    # Freeze first row if requested
    freezePane(wb, sheet = sheet_name, firstRow = freeze_first_row, firstCol = FALSE)
    
    # Bold header
    headerStyle <- createStyle(textDecoration = "bold")
    addStyle(wb, sheet = sheet_name, headerStyle, rows = 1, 
             cols = 1:ncol(ann_ct), gridExpand = TRUE)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  saveWorkbook(wb, file = file.path(output_dir, output_file), overwrite = TRUE)
}

# NTC limma results
ann <- read_rds("../Downloads/limma_perCTex.vivovsin.vivo.rds")
ann <- data.frame(ann)
export_by_celltype(
  df = ann,
  output_dir = here(),
  output_file = "Supplementary_Table1_DEG_cultureEffect.xlsx",
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
  output_dir = here(),
  output_file = "Supplementary_Table2_GSEA_cultureEffect.xlsx",
  sheet_columns = c("pathway", "padj", "NES",  "celltype", "leadingEdge", "db")
)

# KO limma results
ann <- read_rds("../Downloads/limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds")

ann <- ann %>%
  filter(!(coef %in% c("ex.vivo","X.Intercept."))) %>%
  mutate(
    genotype = case_when(
      str_detect(coef, "^ex\\.vivo") ~ str_remove(coef, "^ex\\.vivo"),
      str_detect(coef, "^in\\.vivo") ~ str_remove(coef, "^in\\.vivo"),
      str_detect(coef, "^interaction") ~ str_remove(coef, "^interaction"),
      TRUE ~ NA_character_
    ),
    condition = case_when(
      str_detect(coef, "^ex\\.vivo") ~ "ex.vivo",
      str_detect(coef, "^in\\.vivo") ~ "in.vivo",
      str_detect(coef, "^interaction") ~ "interaction",
      TRUE ~ NA_character_
    ),
    genes = ensg
  )
head(ann)

export_by_celltype(
  df = ann,
  output_dir = here(),
  output_file = "Supplementary_Table3_DEG_interactionEffect.xlsx",
  sheet_columns = c("genotype", "condition", "genes", "logFC", "adj.P.Val", "celltype", "group")
)

# KO enrichment
ann <- fread("../Downloads/fgsea_MSigDB_Hallmark_2020.tsv")
setDT(ann)
ann <- ann[db == "MSigDB_Hallmark_2020"]
ann <- ann[!is.na(NES)]
ann <- data.frame(ann)
ann$condition <- "interaction"

export_by_celltype(
  df = ann,
  output_dir = here(),
  output_file = "Supplementary_Table4_GSEA_interactionEffect.xlsx",
  sheet_columns = c("pathway", "genotype","condition", "padj", "NES", "celltype", "leadingEdge", "db")
)

# culture effect to KO effect
ann <- read_rds("../Downloads/enrichment_to_NTC_genes.rds")
ann <- data.frame(ann)
ann <- ann %>%
  filter(overlap > 5)
export_by_celltype(
  df = ann,
  output_dir = here(),
  output_file = "Supplementary_Table5_overlap_regulatorTargets.xlsx",
  sheet_columns = c("coef", "celltype", "overlap", "p.value", "odds.ratio")
)

# Diff_exp_JAK_STAT
ann <- read_rds("../Downloads/combined_jakstat_diff_exp.rds")
ann <- data.frame(ann)
unique(ann$coef)
ann <- ann %>%
  mutate(
    genotype = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ str_remove(coef, "^ex\\.vivo_genotype"),
      str_detect(coef, "^genotype.*:treatment") ~ str_remove(str_extract(coef, "genotype[^:]+"), "genotype"),
      str_detect(coef, "^genotype") ~ str_remove(coef, "^genotype"),
      str_detect(coef,"treatmentex_vivo") ~ "Wildtype",
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
  output_dir = here(),
  output_file = "Supplementary_Table6_DEG_bulkDataset.xlsx",
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
      str_detect(coef,"treatmentex_vivo") ~ "Wildtype",
      TRUE ~ NA_character_
    ),
    condition = case_when(
      str_detect(coef, "^ex\\.vivo_genotype") ~ "ex.vivo",
      str_detect(coef, ":treatment") ~ "interaction",
      str_detect(coef, "^genotype") ~ "in.vivo",
      TRUE ~ NA_character_
    )
  )
ann <- ann %>%
  filter(condition == "interaction")
export_by_celltype(
  df = ann,
  output_dir = here(),
  output_file = "Supplementary_Table7_GSEA_bulkDataset.xlsx",
  sheet_columns = c("pathway", "genotype", "condition", "padj", "NES", "celltype", "leadingEdge", "db")
)

