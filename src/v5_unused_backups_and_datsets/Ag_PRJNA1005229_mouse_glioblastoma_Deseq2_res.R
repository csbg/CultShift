library(dplyr)
library(readr)
library(stringr)
library(paletteer)
library(ggthemes)
out <- dirout("PRJNA1005229_mouse_glioblastoma_Deseq2_res")
# Set your working directory to where the files are, or provide full path

dir_path <- "/vscratch/aarathy/external_datasets/glioblastoma/"
# List all CSV files in that directory (no subfolders)
files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

# 1. Set directory where CSVs are
dir_path <- "/vscratch/aarathy/external_datasets/glioblastoma/"

# 2. List all CSV files (full paths)
files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

# 3. Keep only noRT files
files_noRT <- files[grepl("_noRT_", files)]

# 4. Extract gene name and condition
file_info <- data.frame(
  filename = files_noRT,
  condition = ifelse(grepl("invitro", files_noRT), "invitro", "preinf"),
  gene = str_extract(basename(files_noRT), "(?<=invitro_|preinf_)[A-Za-z0-9]+(?=_)"),
  stringsAsFactors = FALSE
)

# 5. Initialize list to store merged gene data
merged_data_list <- list()
g <- file_info$gene[1]
# 6. Loop over genes
for(g in unique(file_info$gene)) {
  gene_files <- file_info %>% filter(gene == g)
  dfs <- list()
  for(i in 1:nrow(gene_files)) {
    f <- gene_files$filename[i]
    cond <- gene_files$condition[i]
    
    df <- read.table(f, header = TRUE, sep = " ", quote = "\"", stringsAsFactors = FALSE)
    
    # Numeric indices for all columns except 'feature'
    cols_to_prefix <- which(names(df) != "feature")
    names(df)[cols_to_prefix] <- paste0(cond, "_", names(df)[cols_to_prefix])
    
    dfs[[i]] <- df
  }
  
  merged_gene <- Reduce(function(x, y) full_join(x, y, by = "feature"), dfs)
  merged_gene$gene <- g
  merged_data_list[[g]] <- merged_gene
}

# 7. Combine all genes into one master table
merged_all_genes <- bind_rows(merged_data_list)

# 8. Save as CSV
write_csv(merged_all_genes, "merged_all_genes_noRT.csv")
# 6. Combine all genes into one master table
merged_all_genes <- bind_rows(merged_data_list)

# 7. Save master table
write_csv(merged_all_genes, out("merged_all_genes_noRT.csv"))
###


merged_filtered <- merged_all_genes %>%
  filter(!is.na(preinf_log_fc) & !is.na(invitro_log_fc))





ggplot(merged_filtered, aes(x = preinf_log_fc, y = invitro_log_fc)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  facet_wrap(~ gene, scales = "free") +  # one scatter plot per gene
  labs(
    title = "Scatter plot of in vitro vs pre-infection logFC by gene",
    x = "Pre-infection logFC",
    y = "In vitro logFC"
  ) +
  theme_minimal(base_size = 12)
ggsave(out("scatter_plots.pdf"))
###########
#correlation plots
merged_filtered <- merged_filtered %>%
  mutate(
    invitro_deg = abs(invitro_log_fc) > 1 & invitro_padj < 0.05,
    preinf_deg  = abs(preinf_log_fc)  > 1 & preinf_padj  < 0.05
  )
deg_counts <- merged_filtered %>%
  group_by(gene) %>%
  summarise(
    num_degs = pmax(
      sum(invitro_deg, na.rm = TRUE),
      sum(preinf_deg, na.rm = TRUE)
    ),
    correlation = cor(preinf_log_fc, invitro_log_fc, method = "pearson")
  )



# Get a diverging palette with 30 colors
my_palette <- paletteer_c("grDevices::Purple-Green", 30)
# Example: use it in a barplot or fill scale
ggplot(deg_counts, aes(x = gene, y = correlation, fill = num_degs)) +
  geom_col(color = "black") +
  scale_fill_gradientn(colors = my_palette, name = "num of DEGs") +
  labs(
    x = "Gene",
    y = "Pearson correlation",
    title = "Correlation of in vitro vs pre-infection logFC per gene"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
