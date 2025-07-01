
###############################################

############################
# Filter significant genes
significant_genes <- limmaRes %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

# Classify genes based on log fold changes
significant_genes <- significant_genes %>%
  mutate(
    classification = case_when(
      logFC > 1 & logFC_IRF1KO > 1 ~ "both up",
      logFC < -1 & logFC_IRF1KO < -1 ~ "both down",
      (logFC > 1 & logFC_IRF1KO < -1) | (logFC < -1 & logFC_IRF1KO > 1) ~ "opposite",
      TRUE ~ "other"
    )
  )

# Create the volcano plot
ggplot(significant_genes, aes(x = logFC, y = -log10(adj.P.Val), color = classification)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("both up" = "blue", "both down" = "red", "opposite" = "green", "other" = "black")) +
  labs(
    title = "Volcano Plot of Significant Genes",
    x = "Log Fold Change (Ut_IRF2KO vs Ut_IRF1KO)",
    y = "-log10(Adjusted P-value)"
  ) +
  theme_minimal()
#####################
#
#library(dplyr)

# Filter the significant genes for each condition
significant_IRF2KO <- results %>%
  filter(coef == "Ut_IRF2KO",
         abs(logFC) > 1,
         adj.P.Val < 0.05)

significant_IRF1KO <- results %>%
  filter(coef == "Ut_IRF1KO",
         abs(logFC) > 1,
         adj.P.Val < 0.05)

# Find the overlapping genes
overlap_genes <- inner_join(significant_IRF2KO, significant_IRF1KO, by = "ensg")


# Print the result

# Select relevant columns and classify direction of change
overlap_genes <- overlap_genes %>%
  select(ensg, logFC_IRF2KO = logFC.x, logFC_IRF1KO = logFC.y) %>%
  mutate(direction_IRF2KO = ifelse(logFC_IRF2KO > 1, "up", "down"),
         direction_IRF1KO = ifelse(logFC_IRF1KO > 1, "up", "down"),
         classification = case_when(
           direction_IRF2KO == "up" & direction_IRF1KO == "up" ~ "both up",
           direction_IRF2KO == "down" & direction_IRF1KO == "down" ~ "both down",
           direction_IRF2KO != direction_IRF1KO ~ "opposite"
         ))


# Count the number of genes in each classification category
classification_counts <- overlap_genes %>%
  group_by(classification) %>%
  summarise(count = n())

# Plot as bar chart
ggplot(classification_counts, aes(x = classification, y = count, fill = classification)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Overlap of Significant Genes between Ut_IRF2KO and Ut_IRF1KO",
    x = "Classification",
    y = "Number of Genes"
  ) +
  scale_fill_manual(values = c("both up" = "aquamarine4", 
                               "both down" = "blue4", "opposite" = "darkmagenta")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
base<-"IRF2"
out <- dirout("IRF2")
ggsave(out("Up_down_Opposite.pdf"))
####################################
opposite<-overlap_genes%>%
  filter(classification == "opposite")%>%
  pull(ensg)




opposite <- intersect(opposite, rownames(dataVoom$E))
dat.list<-list()
for(gg in unique(opposite)) {
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
  ggplot(data[data$gene == gene ,], aes(x = sample, y = E)) + 
    geom_boxplot() +
    geom_jitter(aes(colour=sample)) +
    #facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  dir<-dirout(paste0("IRF2"))
  
  ggsave(dir(paste0(gene,".pdf")))
}
colnames(dat.list)
#For each KO where there were significantly up genes, make the plot
for (gg in opposite){
  data<-dat.list%>%filter(gene==gg)%>%filter(treatment_time=="Ut")
 
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
    
  })
}
unique(limmaRes$coef)

head(dat.list)
# Pivot the data to have samples as columns and genes as rows
heatmap_data <- dat.list %>%
  filter(treatment_time == "Ut") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                     "Ut_IRF1_1" ,"Ut_IRF1_2", "Ut_IRF1_3",
                                     "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"))) %>%
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
out
pdf(out("heatmap_Opposite_genes.pdf"),h=16)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "Expression Heatmap",
         fontsize = 8)
dev.off()
##############
#synergy
synergy<-limmaRes%>%filter(coef =="Ut_IRF1KO.IRF2KO")
synergy_genes<- synergy%>%
  filter(
  abs(logFC)>1,adj.P.Val< 0.05)%>%
  pull(ensg)

intersect<-intersect(synergy_genes,opposite)
unique(limmaRes$coef)
