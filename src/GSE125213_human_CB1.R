#Aarathy
###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(readxl)
InDir <- dirout("Ag_top_pathway_genes")
out <- dirout("GSE125213_human_CB1")

uncultured <- read.delim("/media/AGFORTELNY/PROJECTS/TfCf_AG/GSE125213_human_CB/GSE125345_LTHSC_ST_HSC_CMP_GMP_MLP.tsv/GSE125345_LTHSC_ST_HSC_CMP_GMP_MLP.tsv")
uncultured <- as.data.frame(uncultured)
rownames(uncultured) <- make.unique(uncultured$NAME)

gmap <- uncultured[,c("NAME","ENSEMBL_ID")]
uncultured <- read_excel("/media/AGFORTELNY/PROJECTS/TfCf_AG/GSE125213_human_CB/GSE125345_LTHSC_ST_HSC_CMP_GMP_MLP.tsv/GSE125345_countdata.xlsx")
uncul <- as.data.frame(uncultured)
rownames(uncul) <- uncul$ID
uncul$ID <- NULL
# create a named vector: ENSEMBL_ID -> NAME
# create ENSEMBL -> NAME mapping
id2name <- setNames(gmap$NAME, gmap$ENSEMBL_ID)

# get gene names for each row, may contain NAs
gene_names <- id2name[rownames(uncul)]

# remove rows with NA mapping
uncul <- uncul[!is.na(gene_names), ]
gene_names <- gene_names[!is.na(gene_names)]

# drop duplicates (keep first)
uncul <- uncul[!duplicated(gene_names), ]

# set rownames
rownames(uncul) <- gene_names[!duplicated(gene_names)]

# remove temporary column if exists
uncul$gene <- NULL

head(uncul)

boxplot(uncul, las=2)
list <- list.files("/media/AGFORTELNY/PROJECTS/TfCf_AG/GSE125213_human_CB/GSE125213_RAW/")
# list all files
files <- list.files(path = "/media/AGFORTELNY/PROJECTS/TfCf_AG/GSE125213_human_CB/GSE125213_RAW/",
                    pattern = "gene_count.txt.gz$", full.names = TRUE)

# keep only DMSO samples
dmso_files <- grep("DMSO", files, value = TRUE)
read.delim(gzfile(dmso_files[1]), header = TRUE, stringsAsFactors = FALSE)
# read and combine
# read all DMSO files into a list
dmso_list <- lapply(dmso_files, function(f) {
  dat <- read.delim(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
  colnames(dat) <- c("gene", sub("_gene_count.txt.gz", "", basename(f)))
  dat
})

# merge on 'gene'
dmso_counts <- Reduce(function(x, y) merge(x, y, by = "gene"), dmso_list)

# set gene names as rownames
rownames(dmso_counts) <- dmso_counts$gene
dmso_counts <- dmso_counts[,-1]


dmso_counts <- dmso_counts[-c(1:5), ]


# add rownames as a column for merging
dmso_df <- cbind(gene = rownames(dmso_counts), dmso_counts)
uncul_df <- cbind(gene = rownames(uncul), uncul)

# merge by gene
combined_df <- merge(dmso_df, uncul_df, by = "gene")

# set gene names as rownames and remove 'gene' column
rownames(combined_df) <- combined_df$gene
combined_df$gene <- NULL

# final combined matrix
combined_counts <- as.data.frame(combined_df)

#######################
# get sample names
samples <- colnames(combined_counts)

# tissue: DMSO samples are ex vivo, others in vivo
tissue <- ifelse(grepl("Day2_DMSO", samples), "ex.vivo_2d", 
                 ifelse(grepl("Day4_DMSO", samples), "ex.vivo_4d",
                 "in.vivo"))

# celltype: extract prefix before first underscore
# for DMSO samples, it seems CB1/CB2/CB3 info is after first underscore, so we want the part before "_CB"
celltype <- sapply(samples, function(x) {
  if(grepl("DMSO", x)) {
    sub("_CB.*", "", x)   # keeps the part before "_CB"
  } else {
    sub("_CB.*", "", x)   # same for in vivo samples
  }
})

# create metadata dataframe
metadata <- data.frame(
  sample = samples,
  tissue = tissue,
  stringsAsFactors = FALSE
)
metadata$tissue <- factor(metadata$tissue,
                          levels = c("in.vivo","ex.vivo_2d","ex.vivo_4d"))
metadata$celltype <- sapply(1:nrow(metadata), function(i) {
  
  # Take the sample name and tissue for the current row
  sample_name <- metadata$sample[i]
  tissue_type <- metadata$tissue[i]
  
  if (tissue_type == "ex.vivo_2d") {
    # Example: "GSM3565328_CB1_Day2_DMSO"
    # We want to pull out "Day2" (or "Day4", etc.)
    day <- "2d"
    
    # Add "_ex" → becomes "Day2_ex", "Day4_ex"
    celltype <- paste0(day, "_ex")
    
  } 
  if (tissue_type == "ex.vivo_4d") {
    # Example: "GSM3565328_CB1_Day2_DMSO"
    # We want to pull out "Day2" (or "Day4", etc.)
    day <- "4d"
    
    # Add "_ex" → becomes "Day2_ex", "Day4_ex"
    celltype <- paste0(day, "_ex")
    
  }else {
    # For non-ex.vivo, just take everything before "_CB"
    # Example: "Neuron_CB1" → becomes "Neuron"
    celltype <- sub("_CB.*", "", sample_name)
  }
  
  return(celltype)
})
# check result
metadata[, c("sample", "tissue", "celltype")]
design <- model.matrix(~tissue, data=metadata)

#########
d0 <- DGEList(combined_counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]

# voom transformation
dataVoom <- voom(d0, design, plot=TRUE)

#########################

########################
# -----------------------------
# limma fit
# -----------------------------
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

# Extract results
limmaRes <- map_dfr(colnames(coef(limmaFit)), function(coefx) {
  topTable(limmaFit, coef = coefx, number = Inf) %>%
    rownames_to_column("genes") %>%
    filter(coefx != "(Intercept)") %>%
    mutate(coef = coefx,
           group = case_when(
             logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
             logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
             TRUE ~ "n.s"
           ))
})
unique(limmaRes$coef)
#saving D.E table----------------
coef_names <- setdiff(colnames(coef(limmaFit)), "(Intercept)")

limmaRes <- map_dfr(coef_names, function(coefx) {
  
  topTable(limmaFit, coef = coefx, number = Inf, sort.by = "none") %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = coefx,   # 👈 important for your export function
      group = case_when(
        logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
        logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
        TRUE ~ "n.s"
      )
    ) %>%
    dplyr::select(
      gene,
      celltype,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      group
    )
})
export_by_celltype <- function(df, 
                               output_dir, 
                               output_file, 
                               sheet_columns = NULL, 
                               freeze_first_row = TRUE) {
  
  if(!is.null(sheet_columns)){
    df <- df[, intersect(sheet_columns, colnames(df)), drop = FALSE]
  }
  
  df <- as.data.frame(df)
  
  # 👇 KEY CHANGE
  cell_types <- unique(df$celltype)
  split_sheets <- length(cell_types) > 1
  
  wb <- createWorkbook()
  
  if(split_sheets){
    
    for(ct in cell_types){
      ann_ct <- df %>% filter(celltype == ct)
      
      if("adj.P.Val" %in% colnames(ann_ct)){
        ann_ct$adj.P.Val <- format_padj(ann_ct$adj.P.Val)
      }
      ann_ct <- format_numbers(ann_ct)
      
      sheet_name <- substr(gsub("[\\/:*?\\[\\]]", "_", ct), 1, 31)
      
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, ann_ct, rowNames = FALSE)
      
      freezePane(wb, sheet = sheet_name, firstRow = freeze_first_row)
      
      headerStyle <- createStyle(textDecoration = "bold")
      addStyle(wb, sheet = sheet_name, headerStyle,
               rows = 1, cols = 1:ncol(ann_ct), gridExpand = TRUE)
    }
    
  } else {
    # 👇 SINGLE SHEET MODE
    
    ann <- df
    
    if("adj.P.Val" %in% colnames(ann)){
      ann$adj.P.Val <- format_padj(ann$adj.P.Val)
    }
    ann <- format_numbers(ann)
    
    addWorksheet(wb, sheetName = "Results")
    writeData(wb, sheet = "Results", ann, rowNames = FALSE)
    
    freezePane(wb, sheet = "Results", firstRow = freeze_first_row)
    
    headerStyle <- createStyle(textDecoration = "bold")
    addStyle(wb, sheet = "Results", headerStyle,
             rows = 1, cols = 1:ncol(ann), gridExpand = TRUE)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  saveWorkbook(wb, file = file.path(output_dir, output_file), overwrite = TRUE)
}
export_by_celltype(
  df = limmaRes,
  output_dir = out("DE_tables"),
  output_file = "Supplementary_Table5_DE_CB1_limma.xlsx",
  sheet_columns = colnames(limmaRes)
)

###################################
#ranklist
coef_names <- setdiff(colnames(coef(limmaFit)), "(Intercept)")


# Directory for rank lists
dir.create(out("Ranklists"), showWarnings = FALSE)

ranklists <- map(coef_names, function(coefx) {
  
  df <- topTable(
    limmaFit,
    coef = coefx,
    number = Inf,
    sort.by = "none"
  ) %>%
    rownames_to_column("gene")
  
  # Named vector: logFC
  ranks <- df$logFC
  names(ranks) <- df$gene
  
  # Remove NA + duplicates
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  
  # Sort decreasing (important for GSEA)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Save as RDS (best for R)
  saveRDS(
    ranks,
    file = file.path(out("Ranklists"), paste0(coefx, "_logFC_ranklist.rds"))
  )
  
  # Save as TXT (universal)
  write.table(
    data.frame(gene = names(ranks), logFC = ranks),
    file = file.path(out("Ranklists"), paste0(coefx, "_logFC_ranklist.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ranks
})

names(ranklists) <- coef_names
# -----------------------------
# Pivot to long format
# -----------------------------
longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample")%>%
  mutate(tissue_celltype =paste0(tissue,"_",celltype))
longer_dataVoom$organ <- "cord-blood"
longer_dataVoom$author <- "Xie SZ et al., Cell Stem Cell, 2019"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("CB1_NTCs_dataVoom.rds"))

long_data <- as.data.frame(dataVoom$E) %>%
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "sample", values_to = "Expression") %>%
  inner_join(metadata, by = c("sample" = "geo_accession"))

# -----------------------------
# Z-score calculation for gene sets
# -----------------------------
compute_zscore <- function(df, genes) {
  df %>%
    filter(genes %in% genes) %>%
    group_by(genes) %>%
    mutate(
      mean_expr = mean(Expression, na.rm = TRUE),
      sd_expr   = sd(Expression, na.rm = TRUE),
      zscore    = (Expression - mean_expr)/sd_expr
    ) %>%
    ungroup()
}
#####################
ggplot(limmaRes,aes(
  x = logFC,
  y = -log10(adj.P.Val)))+
  geom_point()+
  facet_grid(cols = vars(coef))

#############
ranked_lists <- limmaRes %>%
  
  # Aggregate duplicates: take mean t-statistic per gene per coef
  group_by(coef,genes) %>%
  summarise(logFC = mean(logFC), .groups = "drop") %>%
  
  # Create named vector of ranks per coef
  group_by(coef) %>%
  summarise(
    ranks = list(setNames(logFC, genes)),
    .groups = "drop"
  )


pathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)


fgsea_results <- ranked_lists %>%
  mutate(
    fgsea = purrr::map(ranks, ~ fgsea(pathways = pathways,
                                      stats = .))   # no nperm
  ) %>%
  dplyr::select(coef, fgsea) %>%
  unnest(fgsea)
fgsea_results$dataset <-"CB1"
fgsea_results$organ <-"CB1"
fgsea_results$author <- "Xie SZ et al., Cell Stem Cell, 2019"
fgsea_results %>% write_rds(out("NTC_fgsea_CB1.rds"))
write_rds(fgsea_results, out("FGSEA_CB_2.rds"))

terms <- fgsea_results %>%
  filter(padj < 0.05) %>%
  pull(pathway)%>%
  unique()
fgsea_plot <- fgsea_results %>%
  filter(pathway %in% terms)
#####
# Number of pathways to plot
n_pathways <- length(unique(fgsea_plot$pathway))

# Number of coefficients (if multiple coefs are plotted)
n_coefs <- length(unique(fgsea_plot$coef))

# Define dynamic size (tweak multipliers as needed)
height <- max(4, n_pathways * 0.3)   # 0.3 inches per pathway, min 4 in
width  <- max(6, n_coefs * 2)        # 2 inches per coef, min 6 in

ggplot(fgsea_plot, aes(x = reorder(pathway, NES), y = coef,
                       fill = NES, size = -log10(padj))) +
  geom_point(shape = 21, color = "black") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Pathway",
    y = "Coefficient",
    fill = "NES",
    size = "-log10(padj)",
    title = "FGSEA by Coefficient"
  ) +
  optimized_theme_fig()

#######

##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))


genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

genesets$genes <- toupper(genesets$genes)
#
limmaRes$celltype <- limmaRes$coef
# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)

#
limmaRes_ISGs <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "ISGs")
limmaRes_ISG_core <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISG_core) %>%
  mutate(pathway = "ISGs")
limmaRes_Cholesterol <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "Cholesterol")
limmaRes_mTORC1 <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "mTORC1")


# Define your gene sets (character vectors)
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)

# Combine them into a named list
pathway_list <- list(
  ISGs = ISGs,
  ISG_core = ISG_core,
  Cholesterol = Cholesterol,
  mTORC1 = mTORC1
)

# Iterate over each pathway and gene set, build filtered results, and bind them together
limmaRes_all <- map_dfr(names(pathway_list), function(pathway_name) {
  limmaRes %>%
    select(genes, logFC, adj.P.Val, celltype, group) %>%
    filter(genes %in% pathway_list[[pathway_name]]) %>%
    mutate(pathway = pathway_name)
})

# Optionally combine with genesets if you want one big table
limmaRes_all <- bind_rows(genesets, limmaRes_all)


# Function to generate the dot plot for a given geneset
plot_geneset <- function(data, geneset_name) {
  ggplot(data %>% filter(pathway == geneset_name), 
         aes(x = celltype, y = genes,
             color = pmin(3, pmax(-3, logFC)),
             size = pmin(5, -log10(adj.P.Val)))) +
    geom_point() +
    scale_color_gradient2(
      low = "#4C889C",
      mid = "white",
      high = "#D0154E",
      name = TeX("log_{2}(FC)")
    ) +
    scale_size_continuous(
      range = c(0, 1.8),
      name = TeX("$-\\log_{10}(p_{adj})$")
    ) +
    labs(
      title = paste0("Geneset: ", geneset_name),
      y = "Genes",
      x = "Cell Type"
    ) +
    #facet_grid(rows = vars(pathway), scales = "free", space = "free") +
    optimized_theme_fig()
}

# Make three separate plots
ISG <- plot_geneset(limmaRes_all, "ISGs")
ISG_core <- plot_geneset(limmaRes_all, "ISG_core")
mTORC1 <- plot_geneset(limmaRes_all, "mTORC1")
Cholesterol <- plot_geneset(limmaRes_all, "Cholesterol")

# Save each plot as a PNG
ggsave(out("ISG.png"), ISG, width = 8, height = 6, dpi = 300)
ggsave(out("ISG_core.png"), ISG_core, width = 8, height = 6, dpi = 300)
ggsave(out("mTORC1.png"), mTORC1, width = 8, height = 6, dpi = 300)
ggsave(out("Cholesterol.png"), Cholesterol, width = 8, height = 6, dpi = 300)
###################################
longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample")
longer_dataVoom$tissue <- factor(longer_dataVoom$tissue, 
                                 levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))


prefix <- "CB1_human"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
                                        levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))

# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)




# Filter only genes in your ISGs list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% ISGs)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISG z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(
    title = paste0("Z-score of ISGs across tissues"),
    x = "Tissue",
    y = "Z-score"
  ) +
  facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_ISGs.png"))
########################################
# Define your gene sets (character vectors)
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)

longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample")
longer_dataVoom$tissue <- factor(longer_dataVoom$tissue, 
                                 levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% mTORC1) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% mTORC1)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISG z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
#   labs(
#     title = paste0("Z-score of ISGs across tissues"),
#     x = "Tissue",
#     y = "Z-score"
#   ) +
#   facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
#   optimized_theme_fig() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_mTORC1.png"))

###################
#cholesterol

longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% Cholesterol) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% Cholesterol)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISG z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
#   labs(
#     title = paste0("Z-score of ISGs across tissues"),
#     x = "Tissue",
#     y = "Z-score"
#   ) +
#   facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
#   optimized_theme_fig() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_cholesterol.png"))



