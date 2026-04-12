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
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
basedir <- dirout("Figure13_supplementary")
## InDir4 <- dirout("Figure2_Mye")
limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))

unique(limmaRes$coef)

limmaRes_NTC <- limmaRes %>%
  filter(coef == "tissueex.vivo")
limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
                                                    limmaRes$coef, value = T))
summary_interaction_KO <- limmaRes_int %>%
  # group by interaction term
  mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  mutate(KO = gsub("tissueex.vivo.genotype","",coef.x))%>%
  inner_join(summary_interaction_KO, by = c("KO")) %>%
  filter(signif_count > 10) %>%
  group_by(KO)%>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

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
#order
correlation_results <- correlation_results %>%
  mutate(KO = factor(KO, levels = KO_order))
Correlation_heatmap <- ggplot(correlation_results,
                              aes(x = KO,
                                  y = "Correlation",
                                  fill = cor_abs,
                                  size = pmin(-log10(p_adj), 5))) +
  
  geom_point(shape = 21) +
  
  scale_fill_gradient2(
    low = "purple",
    mid = "white",
    high = "#00A600",
    limits = c(-1, 1),
    name = str_wrap("Correlation", width = 20)
  ) +
  
  scale_size_continuous(
    range = c(0, 2),
    name = TeX("$\\log_{10}\\; (Padj)$")
  ) +
  
  labs(
    title = str_wrap("Correlation of interaction effect to culture effect", width = 80),
    x = "KO",
    y = NULL
  ) +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

Correlation_heatmap

ggsave(basedir("correlation_KO_int_to_ntc.pdf"),
       w= 12, h = 4, units = "cm")
###############

res <- read.delim(paste0(out,"/DEG_ResultsT8.tsv"))

res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")

# Prepare ex.vivo data
res_NTC <- res %>%
  filter(coef == "treatmentex_vivo") %>%
  mutate(genes = gene)%>%
  dplyr::select("genes", "logFC", "cell_type",  "adj.P.Val")
Int_spleen <- res %>% filter(grepl("treatmentex_vivo",res$coef))

Int_spleen$group <- ifelse(Int_spleen$logFC >= 1 & 
                             Int_spleen$adj.P.Val <= 0.05, "up", 
                           ifelse(Int_spleen$logFC <= -1 & 
                                    Int_spleen$adj.P.Val <= 0.05, "down", "n.s"))


# Prepare in.vivo data
Int_spleen <- Int_spleen%>%
  mutate(celltype = gsub("T8", "T-cells", cell_type)) %>%
  filter(celltype == "T-cells")%>%
  filter(coef != "treatmentex_vivo")%>%
  mutate(KO = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(KO = gsub("Interaction_","",KO))%>%
  mutate(genes = gene)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
Int_spleen$dataset <-"Spleen"


summary_int_spleen <- Int_spleen %>%
  # group by interaction term
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")

selected_spleen_KOs <- summary_int_spleen %>%
  filter(signif_count > 10) %>%
  pull(KO) %>%
  unique()


merged_data <-  Int_spleen %>%
  filter(KO %in% selected_spleen_KOs) %>%
  inner_join(res_NTC, by = c("genes"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

correlation_results <- merged_data %>%
  
  group_by(KO)%>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

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
#order

Correlation_heatmap <- ggplot(correlation_results,
                              aes(x = KO,
                                  y = "Correlation",
                                  fill = cor_abs,
                                  size = pmin(-log10(p_adj), 5))) +
  
  geom_point(shape = 21) +
  
  scale_fill_gradient2(
    low = "purple",
    mid = "white",
    high = "#00A600",
    limits = c(-1, 1),
    name = str_wrap("Correlation", width = 20)
  ) +
  
  scale_size_continuous(
    range = c(0, 2),
    limits = c(0,5),
    name = TeX("$\\log_{10}\\; (Padj)$")
  ) +
  
  labs(
    title = str_wrap("Correlation of interaction effect to culture effect", width = 80),
    x = "KO",
    y = NULL
  ) +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

Correlation_heatmap

ggsave(basedir("correlation_KO_int_to_ntc_spleen.pdf"),
       w= 5, h = 2.5, units = "cm")
