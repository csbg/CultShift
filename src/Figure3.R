source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(readr)

#directories ------
base <-"Figure3"
basedir <- dirout("Figure3")
## InDir4 <- dirout("Figure2_Mye")
InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")
#

limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
################################################################################
#gene_based analysis

# Step 1: Filter for significant genes (adj.P.Val < 0.05 and abs(logFC) > 1)
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)  # Only significantly altered genes




data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)
data %>% write_rds(basedir("kos_fig3.rds"))
#Fig3B----------------------------


# grouping gene sets specifically for ISG and cholesterol
genes_fig_1 <- read_rds(InDir_genes("genes_fig1.rds"))
colnames(genes_fig_1) <- c("ensg","pathway")
genes_fig_1$pathway <- gsub("mTORC1_or_Cholesterol","mTORC1/Cholesterol", genes_fig_1$pathway)

filtered_genes <- limmaRes %>%
  inner_join(summary_df,by =c("coef","celltype"))%>%
  filter(Count >10)%>%
  #filter(group != "n.s")%>%
  filter(ensg %in% genes_fig_1$ensg, coef %in% koi) %>%
  group_by(celltype, coef) %>%
  merge(genes_fig_1, by = "ensg") %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)  # Keep only valid KOs for the specific cell type
  
# Recode pathways for better labeling


# fig---------
Fig3B <- ggplot(filtered_genes %>%
                   filter(pathway == "ISG core"), aes(x = coef, y = ensg,
                                                      color = pmin(2, pmax(-2, logFC)) ,
                                                      size = pmin(3, -log10(adj.P.Val))
                   )) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0,1.5),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(title = "Interaction effect of ISG core genes",
      x = "KOs",
    y = "Genes")+
  facet_grid(cols = vars(celltype), rows = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "bottom",
    strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    axis.text.x = element_text(angle = 45))
  

Fig3B
#paper--------------
ggsave(
  filename = basedir("Fig3B.pdf"),
  plot = Fig3B,
  width = 16,
  height = 8.2,
  units = "cm"
)


#
#Fig3C----------------------------
#Fig3Ca---------------------------(Correllation)----------
limmaRes_int <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(coef %in% koi)

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
# Step 4: Plot the correlation results with significance annotations
Fig3Ca <- ggplot(correlation_results, aes(x = coef, y = cor_abs)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  geom_text(aes(label = significance), 
            y =  0.7, 
            color = "black",
            size = 1.5) + 
  labs(
    title = "Correlation of interaction effect to culture effect",
    x = NULL)+
   # fill = "Cell Type")+
  ylab(expression(atop("Correlation of interaction effect", "to culture effects")))+
  
  facet_grid(cols = vars(celltype), scales = "free_x", space = "free_x") +
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        panel.spacing = unit(0.02, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
   
Fig3Ca


#Fig3Cb DEG interaction logFC ---------------------------------------------
ggsave(basedir("Fig3Ca.pdf"),plot=Fig3Ca,
       w=12,h=4, units = "cm")
summary_total <- data %>%
  filter(valid_ko)%>%
  group_by(celltype, coef, genotype, valid_ko) %>%
  summarise(Total_Regulated = sum(Count), .groups = "drop")

Fig3Cb <- ggplot(summary_total,aes(
  x = coef,
  y = log10(Total_Regulated)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  facet_grid(cols = vars(celltype), space = "free", scales = "free") +
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
    #strip.text.x = element_blank(),
    #axis.ticks.x = element_blank())


# Display the plot
Fig3Cb

ggsave(basedir("Fig3Cb.pdf"),plot=Fig3Cb,
       w=12,h=4, units = "cm")

# combine with enrichment to ntc
#Fig3Cc-------------

enrichment_ntc_in.vivo <- read_rds(InDir6("enrichment_to_NTC_genes.rds"))
combined <- data %>%
  dplyr::select(coef,celltype)%>%
  distinct()%>%
  left_join(enrichment_ntc_in.vivo, by = c("coef","celltype"))

# Step 1: Combine the correlation data and enrichment data
combined <- combined %>%
  mutate(significance_en = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))
combined <- combined %>%
  mutate(log2.odds.ratio =log2(odds.ratio))

combined <- combined %>%
  mutate(log2.odds.ratio = case_when(
    coef == "Smc3" ~ 0,
    TRUE ~ odds.ratio
  ))

combined <- combined %>%
  mutate(
    log2.or.filtered = ifelse(overlap > 5, pmin(log2.odds.ratio, 7), NA)  # keep value only if overlap>5
  )

Fig3Cc <- ggplot(combined, aes(x = coef, y = log2.or.filtered)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  # significance labels (for overlap > 5 only, already handled by filtering)
  geom_text(
    data = combined %>% filter(overlap > 5),
    aes(label = significance_en),
    y = 6.5,
    color = "black",
    size = 1.5
  ) +
  
  # add "NA" label for overlap <= 5
  geom_text(
    data = combined %>% filter(overlap <= 5),
    aes(label = "NA"),
    y = 1.5,   # adjust position
    color = "black",
    size = 1.5
  ) +
  
  facet_grid(cols = vars(celltype), scales = "free", space = "free_x") +
  labs(
    x = "KOs",
    y = expression(atop("Overlap of genes with in vivo", "KO effects and culture effects")),
    title = "Overlap of genesets with culture effect and in vivo KO effect"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),  # 👈 center facet labels
    panel.spacing = unit(0.02, "cm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


Fig3Cc

ggsave(basedir("Fig3Cc.pdf"),plot=Fig3Cc,
       w=14.5,h=4, units = "cm")
#Fig3C-----------
Fig3C <-  Fig3Ca / Fig3Cb /  Fig3Cc
ggsave(basedir("Fig3C.pdf"),plot=Fig3C,
       w = 16.5,h=12, units = "cm")

#combined-----------
row1 <- (plot_spacer() | Fig3B) +
  plot_layout(widths = c(1, 1.5))

ggsave(basedir("row1.pdf"), plot = row1, w = 18, h = 8, units = "cm" )
ggsave(basedir("row1_title.pdf"), plot = row1, w = 18, h = 11, units = "cm" )



# Create a truly small blank plot
blank_plot <- ggplot() + theme_void()

# Combine with Fig3C
row2 <- plot_grid(
  blank_plot,
  Fig3C,
  rel_widths = c(1, 3),  # Use relative widths
  nrow = 1
)

#ggsave(basedir("row2_height.pdf"), plot = row2, width = 18, height = 11, units = "cm")


ggsave(basedir("row2.pdf"), row2, w = 18, h = 12, units = "cm" )
#ggsave(basedir("row2_height1.pdf"), row2, w = 18, h = 11, units = "cm" )
##################

