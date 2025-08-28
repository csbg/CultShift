source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")


basedir <- dirout("Figure6_Supplementary")
#Supplementary-------------------

InDir2 <- dirout("Ag_top_filtered_genes")
InDir3 <- dirout("Figure1_Supplementary")
#Sup.Fig5B-----------------
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
#exclude fig1 genes
genes_fig1 <- read_rds(InDir2("genes_fig1.rds"))
# dplyr::select significant genes fron interaction, present in results

supp_fig1_genes <- read_rds(InDir3("supp_fig1_genes.rds"))

genes_ntc <- genes_fig1 %>%
  mutate(geneset = pathway)%>%
  dplyr::select(genes,geneset) %>%
  filter(!(geneset %in% "mTORC1_or_Cholesterol"))%>%
  rbind(supp_fig1_genes)
top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")


top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")%>%
  inner_join(ko_flags, by = c("coef","celltype"))%>%
  filter(valid_ko)%>%
  #filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>% 
  filter(geneset != "ISG core")
# Wrap function for facet labels
wrapped_labeller <- labeller(
  geneset = label_wrap_gen(width = 10) # wrap text after ~8 characters
)

Sup.Fig.6A <- ggplot(top_int, 
                     aes(x = coef, y = ensg,
                         color = pmin(1.5, pmax(-1.5, logFC)),
                         size = pmin(3, -log10(adj.P.Val))
                     )) +
  geom_point() + 
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0, 1.4),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially expressed genesets",
    x = "KOs",
    y = "Genes"
  ) +
  facet_grid(
    rows = vars(geneset), 
    cols = vars(celltype),
    scales = "free", 
    space = "free", 
    labeller = wrapped_labeller
  ) +
  theme_bw() +
  optimized_theme_fig() + 
  theme(
    axis.text = element_text(size = 5),
    strip.text.y = element_text(angle = 90, hjust = 0.5),
    strip.text.x = element_text(angle = 90, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(angle = 0, hjust = 0)  # Rotates legend labels
  )

Sup.Fig.6A



ggsave(basedir("Supplementary_Fig6.pdf"), w=18, h= 18, units = "cm")
