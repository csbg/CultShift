unique(limmaRes$coef)
double_KO_24h<-limmaRes%>%
  filter(
    coef =="doubleKO_24h_IFNb",
    abs(logFC) > 1,
    adj.P.Val < 0.05
  )%>%
  arrange(desc(abs(logFC)))%>%
  head(100)

DKO_24_b <- double_KO_24h %>% pull(ensg)

dat.list<-list()
for(gg in unique(DKO_24_b)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  #avg_gene_data <- gene_data %>%
  #group_by(tissue, celltype) %>%
  #summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
head(dat.list)

##########################################
create_gene_plots_NTC <- function(data, gene) {
  ggplot(data[data$gene == gene ,], aes(x = sample, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter(aes(colour=sample)) +
    #facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  dir<-dirout(paste0("IRF2"))
  
  ggsave(dir(paste0(gene,"_synergy.pdf")))
}
colnames(dat.list)
#For each KO where there were significantly up genes, make the plot
for (gg in DKO_24_b){
  data<-dat.list%>%filter(gene==gg)%>%filter(treatment_time=="IFNb_4h")
  
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
    
  })
}
unique(limmaRes$coef)

head(dat.list)
# Pivot the data to have samples as columns and genes as rows
heatmap_data <- dat.list %>%
  filter(treatment_time == "IFNb_4h") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("IFNb_4h_Rosa_1", "IFNb_4h_Rosa_2", "IFNb_4h_Rosa_3",
                                     "IFNb_4h_IRF2_1", "IFNb_4h_IRF2_2", "IFNb_4h_IRF2_3",
                                     "IFNb_4h_IRF1_1" ,"IFNb_4h_IRF1_2", "IFNb_4h_IRF1_3",
                                     "IFNb_4h_IRF12_1", "IFNb_4h_IRF12_2", "IFNb_4h_IRF12_3"))) %>%
  select(gene, sample1, E) %>%
  spread(sample1, E) %>%
  column_to_rownames(var = "gene")

# Plot heatmap using pheatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "Expression Heatmap",
         fontsize = 8)

pdf(out("heatmap_DKO_24_b_genes.pdf"),h=16)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "Expression Heatmap",
         fontsize = 8)
dev.off()
#
###
#synergy
unique(limmaRes$coef)
# Assuming `limmaRes` contains your results
# Step 1: Extract the interaction terms for 24 hours
interaction_terms <- limmaRes %>%
  filter(coef %in% c(
    "IFNb_24h.IRF1KO.IRF2KO",
    "IFNg_24h.IRF1KO.IRF2KO"
  ))

# Step 2: Filter significant genes
# Assuming an adjusted p-value < 0.05 for significance
significant_genes <- interaction_terms %>%
  filter(adj.P.Val < 0.05)

# Print significant genes
print(significant_genes)

# Step 3: Visualize the results
# Volcano plot
ggplot(significant_genes, aes(x = logFC, y = -log10(adj.P.Val), color = coef)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Synergy Effect at 24h",
       x = "Log Fold Change",
       y = "-log10 Adjusted P-Value") +
  scale_color_manual(values = c("blue", "red"))

ggsave(out("synergy_IFNb_24h.png"))
# Optional: Heatmap of logFC values
logFC_matrix <- significant_genes %>%
  select(ensg, coef, logFC) %>%
  pivot_wider(names_from = coef, values_from = logFC) %>%
  column_to_rownames("ensg") %>%
  as.matrix()
table(is.na(logFC_matrix))
logFC_matrix<-logFC_matrix%>%na.omit()
pdf(out("heatmap_synergy_24h.pdf"))
pheatmap(logFC_matrix, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE)

dev.off()