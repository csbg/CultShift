#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")
#source("src/Ag_enrichR_mouse_genes.R")
library("scales")
library(tidyverse)

basedir <- dirout("Figure5_Supplementary_v8")


limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
KOs <- limmaRes %>%
  pull(coef)%>%
  unique()
coef_sign <- limmaRes_NTC %>% filter(genes %in% KOs) %>%
  filter(group != "n.s")
coef_logFC <- limmaRes_NTC %>% filter(genes %in% KOs)

KO <- koi[1]
ct <- unique(meta$celltype)[1]


dataVoom_NTC_in_ex <- read_rds(InDir_NTC("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(InDir_NTC("NTC_meta.rds"))

example <- c("Cebpa","Brd9","Smc3")
example <- intersect(example, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(example)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = dataVoom_NTC_in_ex$E[gg,]) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
  
  dat.list[[gg]] <- gene_data
}
dat.list <- bind_rows(dat.list,.id="gene")
head(dat.list)
#function
unique(dat.list$celltype)
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  
  # Perform Wilcoxon tests by gene & celltype
  stat_table <- data %>%
    group_by(gene, celltype) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(E ~ tissue)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "n.s.",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "n.s."
      )
    ) %>%
    # Set label y-position based on max expression in each facet
    left_join(
      data %>%
        group_by(gene, celltype) %>%
        summarise(y_pos = max(E, na.rm = TRUE) * 1.05, .groups = "drop"),
      by = c("gene", "celltype")
    )
  
  # Create the plot
  plot <- ggplot(data, aes(x = celltype, y = E, color = tissue, group = tissue)) +
    geom_boxplot(
      fill = NA,
      outlier.shape = NA,
      position = position_dodge(width = 0.8),
      size = 0.2
    ) +
    # geom_jitter(
    #   position = position_jitterdodge(
    #     jitter.width = 0.3,
    #     dodge.width = 0.8
    #   ),
    #   alpha = 0.3,
    #   size = 0.5,
    #   show.legend = FALSE
    # ) +
    geom_text(
      data = stat_table,
      aes(x = celltype, y = y_pos, label = significance),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    facet_grid(
      rows = vars(gene),
      cols = vars(celltype),
      space = "free_x",
      scales = "free",
      labeller = labeller(gene = label_wrap_gen(width = 18))
    ) +
    labs(
      title = "Representative gene expression patterns",
      y = "Scaled Gene Expression",
      x = NULL
    ) +
    scale_color_manual(
      values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
      name = "Experimental model",
      labels = c("ex.vivo" = "Ex vivo", "in.vivo" = "In vivo")
    ) +
    optimized_theme_fig() +
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),    # Light grey major gridlines, subtle but visible
      panel.grid.minor = element_blank()
    )
  
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

# -------------------------------
# Prepare data for plotting
# -------------------------------

# Filter and combine genes
example_genes <- c("Cebpa", "Brd9", "Smc3")
combined_data <- purrr::map_dfr(example_genes, function(gg) {
  df <- dat.list %>%
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP", "Gran.P", "Gran."))
  
  if (nrow(df) == 0) {
    warning(paste("No data available for gene:", gg))
    return(NULL)
  }
  
  df$gene <- gg
  return(df)
})

# Relabel if needed
combined_data$gene <- recode(combined_data$gene,
                             "Cebpa" = "Cebpa",
                             "Brd9" = "Brd9",
                             "Smc3" = "Smc3")

# Set factor levels for consistent ordering
combined_data$celltype <- factor(combined_data$celltype,
                                 levels = c("HSC", "MEP.early", "MkP", 
                                            "GMP", "Gran.P", "Gran.", 
                                            "Mono", "Eo.Ba"))

# Create the plot
Sup.Fig5 <- create_gene_plots_NTC(combined_data, "example_NTC_w.o_jitter", remove_guides = FALSE)

# Save the plot
ggsave(basedir("Sup.Fig5_w.o_jitter.pdf"), Sup.Fig3B, width = 11, height = 7, units = "cm")



