# Define the contrast matrix with reversed order
contrast.mt <- makeContrasts(
  Opposite_Ut = IRF2KO - IRF1KO,
  Opposite_IFNb_4h = IRF2KO + treatment_timeIFNb_4h.IRF2KO - IRF1KO - treatment_timeIFNb_4h.IRF1KO,
  Opposite_IFNb_24h = IRF2KO + treatment_timeIFNb_24h.IRF2KO - IRF1KO - treatment_timeIFNb_24h.IRF1KO,
  Opposite_IFNg_4h = IRF2KO + treatment_timeIFNg_4h.IRF2KO - IRF1KO - treatment_timeIFNg_4h.IRF1KO,
  Opposite_IFNg_24h = IRF2KO + treatment_timeIFNg_24h.IRF2KO - IRF1KO - treatment_timeIFNg_24h.IRF1KO,
  levels = colnames(design)
)

# Fit the contrast
limmaFit.contrast <- contrasts.fit(limmaFit, contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)

# Extract the results for each contrast
for (contrast in colnames(contrast.mt)) {
  print(contrast)
  result <- topTable(limmaFit.contrast, coef=contrast, number=Inf) %>%
    rownames_to_column("ensg") %>%
    mutate(coef=contrast)
  assign(paste0("result_", contrast), result)
}

# Combine the results into a single table if needed
limmaRes_IRF2_vs_IRF1 <- bind_rows(
  result_Opposite_Ut,
  result_Opposite_IFNb_4h,
  result_Opposite_IFNb_24h,
  result_Opposite_IFNg_4h,
  result_Opposite_IFNg_24h
)


################
# Load necessary libraries


###########
# Step 1: Filter genes based on criteria
opposite <- limmaRes_IRF1_vs_IRF2 %>%
  filter(
    abs(logFC) > 1,
    adj.P.Val < 0.05
  ) %>%
  pull(ensg)%>%
  unique()# Assuming ensg is the column with gene IDs

# Example: Assuming you have limma results in limmaRes_IRF1_vs_IRF2
metadata <- metadata %>%
  mutate(
    genotype = case_when(
      IRF1 == "KO" & IRF2 == "KO" ~ "IRF12",
      IRF1 == "WT" & IRF2 == "WT" ~ "WT",
      IRF1 == "KO" & IRF2 == "WT" ~ "IRF1",
      IRF1 == "WT" & IRF2 == "KO" ~ "IRF2",
      TRUE ~ "Unknown"  # Handle any unexpected combinations
    )
  )
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
for (coef in unique(limmaRes_IRF2_vs_IRF1$coef)) {
  
  # Filter significant genes for the current coefficient
  significant_genes <- limmaRes_IRF2_vs_IRF1 %>%
    filter(
      coef == coef,  # Replace 'coef' with the actual coefficient name
      abs(logFC) > 2,
      adj.P.Val < 0.05
    ) %>%
    pull(ensg) %>%
    unique()
  
  # Filter dat_subset for the current coefficient and treatment_time
  dat_subset <- dat.list %>%
    filter(treatment_time == gsub("Opposite_","",coef),
           genotype %in% c("IRF1","IRF2"), # Replace with desired treatment_time
           gene %in% significant_genes)
  
  # Pivot the data to have samples as columns and genes as rows
  heatmap_data <- dat_subset %>%
    select(gene, sample1, E) %>%
    spread(sample1, E) %>%
    column_to_rownames(var = "gene")
  
  # Plot heatmap using pheatmap
  pdf(out(paste("heatmap_", coef, ".pdf", sep="")))
  pheatmap(heatmap_data, 
           cluster_rows = TRUE, 
           cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
           color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
           main = paste("Expression Heatmap -", coef),
           fontsize = 8, 
           show_rownames = FALSE)  # Adjust main title
  dev.off()
}



