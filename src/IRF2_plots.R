base<-dirout("IRF2")
limmaRes<-read_rds(base("limmaRes_all.rds"))

#Allconditions
unique(limmaRes$coef)
logFC_data <- limmaRes %>%
  filter(coef %in% c(
    "IRF1KO_Ut","IRF2KO_Ut","doubleKO_Ut",
    grep("h_IFN",unique(limmaRes$coef),value = T)
    
    
  )) %>%
  select(ensg, coef, logFC) %>%
  pivot_wider(names_from = coef, values_from = logFC)

# Remove rows with any NA values
logFC_data <- logFC_data %>%
  drop_na()

# Calculate pairwise correlations
cor_matrix <- cor(logFC_data[, -1], use = "pairwise.complete.obs")

# Create heatmap
plots<-dirout(paste0("IRF2/Plots"))
pdf(plots("correlation_logFC.pdf"))
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE)
dev.off()

######################################
#Condition:Ut
######################################
condition<-"Ut"
# Separate the data for each coefficient
irf1_data <- limmaRes %>%
  filter(coef == "IRF1KO_Ut") %>%
  select(ensg, logFC_IRF1KO_Ut = logFC, adj.P.Val_IRF1KO_Ut = adj.P.Val)

irf2_data <- limmaRes %>%
  filter(coef == "IRF2KO_Ut") %>%
  select(ensg, logFC_IRF2KO_Ut = logFC, adj.P.Val_IRF2KO_Ut = adj.P.Val)

# Merge the data on the gene identifier
volcano <- inner_join(irf1_data, irf2_data, by = "ensg")

# Inspect the merged data
head(volcano)

# Create the group column based on specified conditions
volcano <- volcano %>%
  mutate(group = case_when(
    logFC_IRF1KO_Ut >= 1 & adj.P.Val_IRF1KO_Ut <= 0.05 & logFC_IRF2KO_Ut >= 1 & adj.P.Val_IRF2KO_Ut <= 0.05 ~ "both_up",
    logFC_IRF1KO_Ut <= -1 & adj.P.Val_IRF1KO_Ut <= 0.05 & logFC_IRF2KO_Ut <= -1 & adj.P.Val_IRF2KO_Ut <= 0.05 ~ "both_down",
    logFC_IRF1KO_Ut >= 1 & adj.P.Val_IRF1KO_Ut <= 0.05 & logFC_IRF2KO_Ut <= -1 & adj.P.Val_IRF2KO_Ut <= 0.05 ~ "IRF1_up_IRF2_down",
    logFC_IRF1KO_Ut <= -1 & adj.P.Val_IRF1KO_Ut <= 0.05 & logFC_IRF2KO_Ut >= 1 & adj.P.Val_IRF2KO_Ut <= 0.05 ~ "IRF1_down_IRF2_up",
    TRUE ~ "others"
  ))
volcano<-as.data.frame(volcano)
# Load ggrepel package
# Define the number of top genes to label
num_top_genes <- 6

# Identify top genes based on adjusted p-values
# Function to select top genes based on logFC within each group

top_genes_IRF1 <- volcano %>%
  filter(group != "others") %>%
  filter(group %in% c("IRF1_down_IRF2_up",
                      "IRF1_up_IRF2_down"
))%>%
  arrange(desc(abs(logFC_IRF1KO_Ut))) %>%
  head(num_top_genes)
top_genes_IRF2 <-volcano %>%
  filter(group != "others") %>% 
  filter(group %in% c("IRF1_down_IRF2_up","IRF1_up_IRF2_down")
         )%>%
  arrange(desc(abs(logFC_IRF2KO_Ut))) %>%
  head(num_top_genes)

# Inspect the data after adding the group column

#################################################
# Plot the data
# Function to calculate label positions and endpoints for lines
calculate_label_ends <- function(data, x_offset, y_offset) {
  data <- mutate(data,
                 label_x = logFC_IRF1KO_Ut + x_offset,
                 label_y = logFC_IRF2KO_Ut + y_offset,
                 line_x = ifelse(logFC_IRF1KO_Ut < label_x, logFC_IRF1KO_Ut + 0.2, logFC_IRF1KO_Ut - 0.2),
                 line_y = ifelse(logFC_IRF2KO_Ut < label_y, logFC_IRF2KO_Ut + 0.2, logFC_IRF2KO_Ut - 0.2))
  return(data)
}

# Calculate label end points for each group
top_genes_IRF1 <- calculate_label_ends(top_genes_IRF1, x_offset = 0.5, y_offset = 0.5)
top_genes_IRF2 <- calculate_label_ends(top_genes_IRF2, x_offset = 0.5, y_offset = -0.5)

# Plot the data
ggplot() +
  # Hexbin plot for the "others" group
  stat_bin_hex(data = filter(volcano, group == "others"), 
               aes(x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut, fill = ..count..), 
               bins = 50, color = NA, alpha = 0.7) +
  scale_fill_gradient(low = "lightgrey",  high= "steelblue", limits = c(1, 5000), name = "Gene Count") +
  
  # Overlay points for the main groups
  geom_point(data = filter(volcano, group != "others"), 
             aes(x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut, color = group), 
             alpha = 0.9, size = 2.5) +
  
  # Manually setting colors for groups
  scale_color_manual(values = c(
    "both_up" = "#4C889C", 
    "both_down" = "#4C889C", 
    "IRF1_up_IRF2_down" = "#D0154E", 
    "IRF1_down_IRF2_up" = "#D0154E"
  ), name = "Group") +
  
  # Add labels for top genes with lines using geom_text_repel
  # Add labels for top genes with lines using geom_text_repel
  geom_text_repel(data = top_genes_IRF1,
                  aes(label = ensg, x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut), 
                  size = 3, 
                  box.padding = 0.3, 
                  point.padding = 0.3, 
                  segment.color = 'grey50',
                  nudge_x = -0.9, nudge_y = -0.7,  # Adjust label position
                  force = 5,  # Increase force for stronger repelling
                  min.segment.length = unit(0.2, "cm")) +
  
  geom_text_repel(data = top_genes_IRF2,
                  aes(label = ensg, x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut), 
                  size = 3, 
                  box.padding = 0.3, 
                  point.padding = 0.3, 
                  segment.color = 'grey50',
                  nudge_x = 0.9, nudge_y = -0.7,
                  force = 5,
                  min.segment.length = unit(0.2, "cm")) +
  
  # Add lines from points to labels
  # geom_segment(data = top_genes_IRF1,
  #              aes(x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut, xend = line_x, yend = line_y),
  #              color = "black", alpha = 1, size = 0.5, lineend = "round") +
  # 
  # geom_segment(data = top_genes_IRF2,
  #              aes(x = logFC_IRF1KO_Ut, y = logFC_IRF2KO_Ut, xend = line_x, yend = line_y),
  #              color = "black", alpha = 1, size = 0.5, lineend = "round") +
  
  labs(title = "IRF1KO vs IRF2KO",
       x = "logFC IRF1KO_Ut",
       y = "logFC IRF2KO_Ut") +
  
  theme_minimal() +
  
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

plots<-dirout(paste0("IRF2/",condition))
ggsave(plots("Volcano_",condition,".pdf"))

#################################################
#Antagonistic_based_on_logFC
################################################
gene_set<-"Antagonism"
Antagonistic_Ut <- volcano%>%
  filter(group %in% c("IRF1_down_IRF2_up","IRF1_up_IRF2_down"))%>%
  pull(ensg)%>%
  unique()

metadata <- metadata %>%
  mutate(
    genotype = case_when(
      IRF1 == "KO" & IRF2 == "KO" ~ "IRF12KO",
      IRF1 == "WT" & IRF2 == "WT" ~ "WT",
      IRF1 == "KO" & IRF2 == "WT" ~ "IRF1KO",
      IRF1 == "WT" & IRF2 == "KO" ~ "IRF2KO",
      TRUE ~ "Unknown"  # Handle any unexpected combinations
    )
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "IRF1KO", "IRF2KO", "IRF12KO"))
  )
dat.list<-list()

for(gg in unique(Antagonistic_Ut)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = dataVoom$E[gg,]) %>%
    mutate(scaled_E=scale(E))%>%
    rownames_to_column("sample1") %>%
    filter(genotype %in% c("WT","IRF2KO","IRF1KO","IRF12KO"),treatment_time=="Ut")%>%
    remove_rownames()
  

  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")

##########################################
create_gene_plots_NTC <- function(data, gene,gene_set) {
  ggplot(data[data$gene == gene ,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter(aes(color=genotype)) +
    #facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  dir<-dirout(paste0("IRF2/",gene_set,"/Genes"))
  
  ggsave(dir(paste0(gene,".pdf")))
}

#For each KO where there were significantly up genes, make the plot
for (gg in Antagonistic_Ut){
  data<-dat.list%>%filter(gene==gg)%>%filter(treatment_time=="Ut")
  
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene,"Antagonism")
    
  })
}

# Pivot the data to have samples as columns and genes as rows
heatmap_data <- dat.list %>%
  filter(treatment_time == "Ut") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                     "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                     "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"
                                     ))) %>%
  select(gene, sample1, scaled_E) %>%
  spread(sample1,scaled_E) %>%
  column_to_rownames(var = "gene")

# Plot heatmap using pheatmap

dir<-dirout(paste0("IRF2/",gene_set,"/heatmap"))
pdf(dir("heatmap_Antagonistic_Ut_genes.pdf"),h=length(Antagonistic_Ut)*0.15)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "IRF1-IRF2-antagonism",
         fontsize = 8)
dev.off()
###############################################################################

gene_set<-"Antagonism"
# Ensure coef is in the specified order
Antagonistic <- limmaRes %>%
  filter(coef %in% c("IRF1KO_Ut", "IRF2KO_Ut", "doubleKO_Ut", "IRF2_vs_IRF1_Ut", "IRF1KO.IRF2KO_Ut")) %>%
  mutate(coef = factor(coef,
                       levels = c("IRF1KO_Ut", "IRF2KO_Ut", "doubleKO_Ut", "IRF2_vs_IRF1_Ut", "IRF1KO.IRF2KO_Ut")))%>%
  filter(
    ensg %in% Antagonistic_Ut)

# Custom x-axis labels
custom_labels <- c("IRF1KO_Ut" = "IRF1KO vs Wt", 
                   "IRF2KO_Ut" = "IRF2KO vs Wt", 
                   "doubleKO_Ut" = "DoubleKO vs Wt", 
                   "IRF2_vs_IRF1_Ut" = "IRF2KO vs IRF1KO", 
                   "IRF1KO.IRF2KO_Ut" = "IRF1-IRF2 interaction ")
# Plot
ggplot(Antagonistic, aes(x = coef, y = ensg, size = -log10(adj.P.Val), color = logFC)) +
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_color_gradient2(high = "red", low = "blue") +
  scale_x_discrete(labels = custom_labels) +
  theme(axis.text = element_text(size = 15)) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  #scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
  ggtitle(paste0(gene_set)) 
  #ylab("Antagonistic genes")
ggsave(dir("LogFC_Ut_genes.pdf"),h=length(Antagonistic_Ut)*0.2)


#########################
#Differences between IRF1KO and IRF2KO
##########################
gene_set<-"IRF2_IRF1_difference"
Difference_Ut<-limmaRes%>%
  filter(coef == "IRF2_vs_IRF1_Ut",
         abs(logFC)>1, 
         adj.P.Val< 0.05)%>%
  pull(ensg)
dat.list<-list()
for(gg in unique(Difference_Ut)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = dataVoom$E[gg,]) %>%
    mutate(scaled_E=scale(E))%>%
    rownames_to_column("sample1") %>%
    filter(genotype %in% c("WT","IRF2KO","IRF1KO","IRF12KO"),treatment_time=="Ut")%>%
    remove_rownames()
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}  
dat.list<-bind_rows(dat.list,.id="gene")


heatmap_data <- dat.list %>%
  filter(treatment_time == "Ut") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                     "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                     "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"
                          ))) %>%
  select(gene, sample1, scaled_E) %>%
  spread(sample1,scaled_E) %>%
  column_to_rownames(var = "gene")

# Plot heatmap using pheatmap

dir<-dirout(paste0("IRF2/",gene_set,"/heatmap/"))
pdf(dir(paste0("heatmap.pdf")))
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "IRF1-IRF2-difference",
         fontsize = 8,
         show_rownames = F)
dev.off()


###############################################
#fgsea
#############################################

enr.terms <- enrichrGetGenesets(ENRICHR.DBS)

# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})

##############################################
set.seed(123)  # For reproducibility
gsea.res <- data.table()
# Add a small amount of random noise to logFC to break ties
unique(limmaRes$coef)
dbx <- "MSigDB_Hallmark_2020"
IRF2 <- limmaRes %>%filter(coef =="IRF2KO_Ut")%>%
  mutate(logFC_no_tie = logFC + rnorm(n(), mean = 0, sd = 1e-5))


# Now using `with` correctly on the subsetted data
stats <- with(IRF2, setNames(logFC_no_tie, nm = ensg))

# Check for missing values in stats
if (any(is.na(stats))) {
  next  # Skip this iteration if there are missing values in stats
}

# Perform fgsea analysis
fgsea_output <- fgsea(
  pathways = enr.terms[["MSigDB_Hallmark_2020"]],
  stats = stats
  #minSize = 15,   # Example additional arguments, adjust as necessary
  #maxSize = 500,  # Example additional arguments, adjust as necessary
  #nperm = 1000    # Example additional arguments, adjust as necessary
)

gsea.res <- data.table()  # Initialize the result table

gsea.res<- fgsea_output
# cleanup / export limmaRes
gsea.res[is.nan(NES), NES = 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1.10], collapse = ","))

#
pDT <- gsea.res

## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5)]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5)]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]
#pDT <- hierarch.ordering(pDT, "pathway", "celltype", "NES", TRUE)
#pDT <- hierarch.ordering(pDT, "pathway", "NES",T)
ggplot(pDT, aes(x=size, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
  #geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red") +
  geom_point(data=pDT[padj < 0.05]) +
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  theme_bw(12) +
  xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(x = "No:of genes")+
  
  theme(axis.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#theme(strip.text.y=element_text(angle=0))
dat<-dirout(paste0("IRF2/","FGSEA/",dbx))
ggsave(dat("GSEA_plot_","MSigDb_IRF2",".pdf"),h=5,w=7)


# Function to retrieve leading edge genes
getLeadingEdgeGenes <- function(pathway_name) {
  leading_edge <- leadingEdge(gsea.res, pathway_name)
  return(leading_edge$gene)
}

# Apply function to each pathway in pDT_significant
leading_edge_genes <- lapply(pDT$pathway, getLeadingEdgeGenes)
extract_unique_genes <- function(fgsea,
                                 terms) {
  unique_genes <- fgsea %>%
    filter(pathway %in% terms) %>%
    pull(leadingEdge) %>%
    #unique() %>%
    #strsplit(",") %>%
    unlist() %>%
    #toupper() %>%
    # unique()
    return(unique_genes)
}
unique_pathways<-unique(pDT$pathway)
mapped_genes <- map(unique_pathways, ~ extract_unique_genes(pDT, .))
mapped_genes <- setNames(mapped_genes, unique_pathways)
#############################################################
TGF_genes<-c(mapped_genes$`Interferon Alpha Response`,mapped_genes$`Interferon Gamma Response`)
TGF_genes<-unique(TGF_genes)
TGF_genes<-c(mapped_genes$`TGF-beta Signaling`)

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

for(gg in unique(TGF_genes)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    filter(genotype %in% c("WT","IRF2"),treatment_time=="Ut")%>%
    
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
  
  ggsave(dir(paste0(gene,"_TGF.pdf")))
}
colnames(dat.list)
#For each KO where there were significantly up genes, make the plot
for (gg in TGF_genes){
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
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3"))) %>%
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
pdf(out("heatmap_TGF_genes_genes.pdf"))
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "Expression Heatmap",
         fontsize = 8)
dev.off()
##################################################################
#
unique(limmaRes$coef)
limmaRes_plot<-limmaRes%>%
  filter(
    coef %in% c("IRF2KO_Ut","IRF1KO_Ut","IRF2_vs_IRF1_Ut" ,"IRF1KO.IRF2KO_Ut","doubleKO_Ut")
  )
coef_order <- c("IRF2KO_Ut","IRF1KO_Ut","IRF2_vs_IRF1_Ut" ,"IRF1KO.IRF2KO_Ut","doubleKO_Ut")
top_genes<-limmaRes_plot%>%
  filter(
    abs(logFC) > 1,
    adj.P.Val < 0.05
  ) %>%
  arrange(coef, desc(logFC)) %>%
  group_by(coef) %>%
  slice_head(n = 10) %>%
  ungroup()%>%
  pull(ensg)%>%
  unique()
top_genes <- limmaRes_plot %>%
  filter(ensg %in% top_genes)%>%arrange(factor(coef, levels = coef_order), desc(logFC))
top_genes$coef<-factor(top_genes$coef, levels = coef_order)

ggplot(top_genes, aes(x = reorder(ensg, logFC), y = logFC)) +
  geom_point +
  facet_grid(cols=vars(coef), scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 15)) +
  theme(axis.text.y = scaled_Element_text(size = 15))+
  labs(x = "Genes", y = "logFC", title = "Top 10 Genes by logFC for each coef") +
  coord_flip()


################################################################################
ggplot(top_genes, aes(x = coef, y = ensg, size = -log10(adj.P.Val), color = logFC)) +
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_color_gradient2(high = "red", low = "blue") +
  scale_x_discrete(labels = custom_labels) +
  theme(axis.text = element_text(size = 15)) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  #scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
  ggtitle(paste0("Ut-Top 10 genes")) 
#ylab("Antagonistic genes")
ggsave(plots("top_10_genes_Ut.pdf"),h=length(unique(top_genes$ensg))*0.22)
################################################################################
IRF1_or_IRF2<-limmaRes%>%
  filter(coef %in% c("IRF1KO_Ut","IRF2KO_Ut"),
         abs(logFC)>1,
                adj.P.Val<0.05)%>%
  pull(ensg)%>%
  unique()

gene_set <- "IRF1_or_IRF2" 
dat.list<-list()
for(gg in unique(IRF1_or_IRF2)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = dataVoom$E[gg,]) %>%
    mutate(scaled_E=scale(E))%>%
    rownames_to_column("sample1") %>%
    filter(genotype %in% c("WT","IRF2KO","IRF1KO","IRF12KO"),treatment_time=="Ut")%>%
    remove_rownames()
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}  
dat.list<-bind_rows(dat.list,.id="gene")


heatmap_data <- dat.list %>%
  filter(treatment_time == "Ut") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                     "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                     "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"
                          ))) %>%
  select(gene, sample1, scaled_E) %>%
  spread(sample1,scaled_E) %>%
  column_to_rownames(var = "gene")
dir<-dirout(paste0("IRF2/",gene_set))
pdf(dir("heatmap_Ut_genes.pdf"))
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "IRF1_or_IRF2_regulated",
         fontsize = 8,
         show_rownames = F)
dev.off()
#################################
#########################
#Interaction effect 
########################
gene_set<-"Interaction"
Interaction_IRF1_IRF2<-limmaRes%>%
  filter(coef %in% c("IRF1KO.IRF2KO_Ut"),
         abs(logFC)>1,
         adj.P.Val<0.05)%>%
  pull(ensg)%>%
  unique()


dat.list<-list()
for(gg in unique(Interaction_IRF1_IRF2)) {
  # Subset the metadata and E values for the current gene
  gene_data <- metadata %>%
    mutate(E = dataVoom$E[gg,]) %>%
    mutate(scaled_E=scale(E))%>%
    rownames_to_column("sample1") %>%
    filter(genotype %in% c("WT","IRF2KO","IRF1KO","IRF12KO"),treatment_time=="Ut")%>%
    remove_rownames()
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}  
dat.list<-bind_rows(dat.list,.id="gene")
#genes
create_gene_plots_NTC <- function(data, gene,gene_set) {
  ggplot(data[data$gene == gene ,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter(aes(color=genotype)) +
    #facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  dir<-dirout(paste0("IRF2/",gene_set,"/Genes"))
  
  ggsave(dir(paste0(gene,".pdf")))
}

#For each KO where there were significantly up genes, make the plot
gene_set
for (gg in Interaction_IRF1_IRF2){
  data<-dat.list%>%filter(gene==gg)%>%filter(treatment_time=="Ut")
  
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene,gene_set)
    
  })
}

#heatmap
heatmap_data <- dat.list %>%
  filter(treatment_time == "Ut") %>%
  mutate(sample1 = factor(sample1, 
                          levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                     "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                     "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                     "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"
                          ))) %>%
  select(gene, sample1, scaled_E) %>%
  spread(sample1,scaled_E) %>%
  column_to_rownames(var = "gene")
dir<-dirout(paste0("IRF2/",gene_set))
pdf(dir("heatmap_Ut_genes.pdf"),h=length(Interaction_IRF1_IRF2)*0.1)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "IRF1_IRF2_interaction",
         fontsize = 8,
         show_rownames = T)
dev.off()

##################
#numbers
# Filter the significant genes for each condition
significant_IRF2KO <- limmaRes %>%
  filter(coef == "IRF2KO_Ut",
         abs(logFC) > 1,
         adj.P.Val < 0.05)%>%nrow()

significant_IRF1KO <- limmaRes %>%
  filter(coef == "IRF1KO_Ut",
         abs(logFC) > 1,
         adj.P.Val < 0.05)%>%nrow()

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

tab <- dirout("IRF2/tables")
write.table(overlap_genes,tab("overlap_IRF1_IRF2.tsv"))
# Count the number of genes in each classification category
classification_counts <- overlap_genes %>%
  group_by(classification) %>%
  summarise(count = n())
limmaRes %>%
  filter(coef == "IRF1KO_Ut",
         logFC > 1,
         adj.P.Val < 0.05)%>%nrow()
limmaRes %>%
  filter(coef == "IRF1KO_Ut",
         logFC < -1,
         adj.P.Val < 0.05)%>%nrow()
limmaRes %>%
  filter(coef == "IRF2KO_Ut",
         logFC > 1,
         adj.P.Val < 0.05)%>%nrow()
limmaRes %>%
  filter(coef == "IRF2KO_Ut",
         logFC < -1,
         adj.P.Val < 0.05)%>%nrow()
classification_counts <- classification_counts %>%
  add_row(classification = "IRF1KO_up", count = 534) %>%
  add_row(classification = "IRF2KO_up", count = 168) %>%
  add_row(classification = "IRF1KO_down", count = 514) %>%
  add_row(classification = "IRF2KO_down", count = 121)


# Plot as bar chart
ggplot(classification_counts, aes(x = classification, y = count)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Overlap of Significant Genes between IRF2KO_Ut and IRF1KO_Ut",
    x = "Classification",
    y = "Number of Genes"
  ) +
  #scale_fill_manual(values = c("both up" = "aquamarine4", 
                              # "both down" = "blue4", "opposite" = "darkmagenta")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

plots<-dirout(paste0("IRF2/Plots"))
ggsave(plots("Up_down_Opposite.pdf"))
####################################################
results<-limmaRes%>%
  filter(
    coef %in% c("IRF2KO_Ut","IRF1KO_Ut","IRF2_vs_IRF1_Ut" ,"IRF1KO.IRF2KO_Ut","doubleKO_Ut")
  )
write.table(results,tab("Results_differential_expression.tsv"))
