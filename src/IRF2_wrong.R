
# Load necessary libraries
library(dplyr)
library(limma)
library(pheatmap)
library(here)
library(edgeR)
############################################

# Define the function
get_top_genes <- function(data, coefs_of_interest=unique(data$coef), top_n) {
  data %>%
    filter(coef %in% coefs_of_interest) %>%
    group_by(coef) %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

#
# Load the count data
countdata <- read.csv("/Users/aarathyrg/Downloads/raw_IRF2_gene.csv", header = TRUE)
rownames(countdata) <- countdata$X
countdata$X <- NULL
countdata <- as.matrix(countdata)

# Extract sample names
samples <- colnames(countdata)
treatments <- c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
# Define metadata with correct assignment of treatment_time and IRF genotype
metadata <- data.frame(
  sample = samples,
  treatment_time = factor(
    ifelse(
      grepl("IFNb_4h", samples), "IFNb_4h",
      ifelse(
        grepl("IFNb_24h", samples), "IFNb_24h",
        ifelse(
          grepl("IFNg_4h", samples), "IFNg_4h",
          ifelse(
            grepl("IFNg_24h", samples), "IFNg_24h",
            "Ut"  # Default for untreated samples
          )
        )
      )
    ),
    levels = c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
  ),
  IRF1 = factor(
    ifelse(grepl("IRF1", samples) | grepl("IRF12", samples), "KO", "WT"),
    levels = c("WT", "KO")
  ),
  IRF2 = factor(
    ifelse(grepl("IRF2", samples) | grepl("IRF12", samples), "KO", "WT"),
    levels = c("WT", "KO")
  )
)

# Ensure Rosa is set as wildtype
metadata$IRF1[grepl("Rosa", metadata$sample)] <- "WT"
metadata$IRF2[grepl("Rosa", metadata$sample)] <- "WT"

# Set row names of metadata to sample names
rownames(metadata) <- metadata$sample

# Ensure that the columns in countdata match the row names in metadata
countdata <- countdata[, rownames(metadata)]
stopifnot(all(colnames(countdata) == rownames(metadata)))
#############################################################
# Specify the design formula without an intercept
#############################################################
#edgeR norm
# Create the design matrix
design <- model.matrix(~treatment_time * IRF1 * IRF2, data = metadata)
# Filtering
d0 <- DGEList(countdata)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dir<-"Olga_Limma_Results/"
pdf(paste0(dir,"voom.png"))
dataVoom<-voom(d, design, plot = TRUE)
dev.off()

############################################################
#fit model
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

limmaRes <- list() # start an empty list
for(coefx in colnames(coef(limmaFit))){ # run a loop for each coef
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coef in a list.
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) %>%
    rownames_to_column("ensg")
}
limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coef in the column "coef"
################################################
#conrtrasts
# Define treatment conditions
treatments <- c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")

# Initialize the contrast matrices with zeros
contrast_matrix_IRFDKO_SKOs <- matrix(0, nrow = length(colnames(coef(limmaFit))),
                                      ncol = length(treatments))
contrast_matrix_IRF2KO_vs_IRF1KO <- matrix(0, nrow = length(colnames(coef(limmaFit))),
                                           ncol = length(treatments))
rownames(contrast_matrix_IRFDKO_SKOs) <- colnames(coef(limmaFit))
rownames(contrast_matrix_IRF2KO_vs_IRF1KO) <- colnames(coef(limmaFit))
colnames(contrast_matrix_IRFDKO_SKOs) <- paste0(treatments, "_IRFDKO-SKOs")
colnames(contrast_matrix_IRF2KO_vs_IRF1KO) <- paste0(treatments, "_IRF2KO_vs_IRF1KO")

# Define the contrasts for each treatment condition
for (treatment in treatments) {
  if (treatment == "Ut") {
    contrast_matrix_IRFDKO_SKOs[grep("^IRF1KO:IRF2KO$", rownames(contrast_matrix_IRFDKO_SKOs)),
                                paste0(treatment, "_IRFDKO-SKOs")] <- 1
    contrast_matrix_IRFDKO_SKOs[grep("^IRF1KO$", rownames(contrast_matrix_IRFDKO_SKOs)),
                                paste0(treatment, "_IRFDKO-SKOs")] <- -1
    contrast_matrix_IRFDKO_SKOs[grep("^IRF2KO$", rownames(contrast_matrix_IRFDKO_SKOs)),
                                paste0(treatment, "_IRFDKO-SKOs")] <- -1
    
    contrast_matrix_IRF2KO_vs_IRF1KO[grep("^IRF2KO$",rownames(contrast_matrix_IRF2KO_vs_IRF1KO)),
                                     paste0(treatment, "_IRF2KO_vs_IRF1KO")] <- 1
    contrast_matrix_IRF2KO_vs_IRF1KO[grep("^IRF1KO$", rownames(contrast_matrix_IRF2KO_vs_IRF1KO)),
                                     paste0(treatment, "_IRF2KO_vs_IRF1KO")] <- -1
  } else {
    contrast_matrix_IRFDKO_SKOs[grep(
      paste0(treatment, ":IRF1KO:IRF2KO"), rownames(contrast_matrix_IRFDKO_SKOs)
    ), paste0(treatment, "_IRFDKO-SKOs")] <- 1
    contrast_matrix_IRFDKO_SKOs[grep(
      paste0(treatment, ":IRF1KO$"), rownames(contrast_matrix_IRFDKO_SKOs)
    ), paste0(treatment, "_IRFDKO-SKOs")] <- -1
    contrast_matrix_IRFDKO_SKOs[grep(
      paste0(treatment, ":IRF2KO$"), rownames(contrast_matrix_IRFDKO_SKOs)
    ), paste0(treatment, "_IRFDKO-SKOs")] <- -1
    
    contrast_matrix_IRF2KO_vs_IRF1KO[grep(
      paste0(treatment, ":IRF2KO$"), rownames(contrast_matrix_IRF2KO_vs_IRF1KO)
    ), paste0(treatment, "_IRF2KO_vs_IRF1KO")] <- 1
    contrast_matrix_IRF2KO_vs_IRF1KO[grep(
      paste0(treatment, ":IRF1KO$"), rownames(contrast_matrix_IRF2KO_vs_IRF1KO)
    ), paste0(treatment, "_IRF2KO_vs_IRF1KO")] <- -1
  }
}


# Fit the contrasts to the limma model
limmaFit_contrast_IRFDKO_SKOs <- contrasts.fit(limmaFit, contrast_matrix_IRFDKO_SKOs)
limmaFit_contrast_IRFDKO_SKOs <- eBayes(limmaFit_contrast_IRFDKO_SKOs)

limmaFit_contrast_IRF2KO_vs_IRF1KO <- contrasts.fit(limmaFit, contrast_matrix_IRF2KO_vs_IRF1KO)
limmaFit_contrast_IRF2KO_vs_IRF1KO <- eBayes(limmaFit_contrast_IRF2KO_vs_IRF1KO)


limmaRes.contrast_IRFDKO_SKOs_list <- lapply(colnames(contrast_matrix_IRFDKO_SKOs),
                                             function(coef_name) {
                                               topTable(limmaFit_contrast_IRFDKO_SKOs, coef = coef_name, number = Inf) %>%
                                                 rownames_to_column("ensg") %>%
                                                 mutate(coef = coef_name)
                                             })
# Combine all contrast results into a single data frame
limmaRes.contrast_IRFDKO_SKOs <- do.call(rbind, limmaRes.contrast_IRFDKO_SKOs_list)

limmaRes.contrast_IRF2KO_vs_IRF1KO_list <- lapply(colnames(contrast_matrix_IRF2KO_vs_IRF1KO),
                                                  function(coef_name) {
                                                    topTable(limmaFit_contrast_IRF2KO_vs_IRF1KO,
                                                             coef = coef_name, number = Inf) %>%
                                                      rownames_to_column("ensg") %>%
                                                      mutate(coef = coef_name)
                                                  })
# Combine all contrast results into a single data frame
limmaRes.contrast_IRF2KO_vs_IRF1KO <- do.call(rbind, limmaRes.contrast_IRF2KO_vs_IRF1KO_list)


# Add them to the full table (limmaRes)
limmaRes <- rbind(limmaRes.contrast_IRF2KO_vs_IRF1KO,
                  limmaRes.contrast_IRFDKO_SKOs, limmaRes)

################################################
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept

############################
limmaRes$coef<-gsub("treatment_time","",limmaRes$coef)

limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    coef %in% c("IRF1KO", "IRF2KO", "IRF1KO:IRF2KO") ~ paste0("Ut_", coef),
    TRUE ~ coef
  ))
############################
#adding treatment and genotype
# Extract treatment and genotype from coef using regex
limmaRes <- limmaRes %>%
  mutate(
    treatment = ifelse(grepl("IFNb_4h", coef), "IFNb_4h",
                       ifelse(grepl("IFNg_4h", coef), "IFNg_4h",
                              ifelse(grepl("IFNb_24h", coef), "IFNb_24h",
                                     ifelse(grepl("IFNg_24h", coef), "IFNg_24h",
                                            "Ut")))),
    genotype = ifelse(grepl("IRF1KO:IRF2KO", coef), "IRF12KO",
                      ifelse(grepl("IRF2KO_vs_IRF1KO", coef), "IRF2KO_vs_IRF1KO",
                             ifelse(grepl("IRF1KO", coef), "IRF1KO",
                                    ifelse(grepl("IRF2KO", coef), "IRF2KO",
                                           ifelse(grepl("IRFDKO-SKOs", coef), "IRFDKO-SKOs",
                                                  
                                                  "WT")
                                    )
                             )
                      )
    )
  )%>%
  mutate(color = ifelse(adj.P.Val < 0.05, "Significant", "N.S"))



###############################################################################
unique(limmaRes$coef)
genes<-limmaRes%>%
  filter(coef == "Ut_IRF1KO:IRF2KO")%>%
  filter(adj.P.Val < 0.05, abs(logFC)>1 )%>%
  arrange(desc(abs(logFC)))%>%
  head(50)%>%
  pull(ensg)%>%
  unique()
#top 20 genes regulate by IRF1KO/IRF2 KO at untreated
IRF1_IRF2_DKO<-limmaRes%>%
  filter(ensg %in% genes,
         coef %in% c("Ut_IRF1KO","Ut_IRF2KO","Ut_IRF1KO:IRF2KO"),
         adj.P.Val< 0.05)
IRF1_IRF2_DKO$adj.P.Val<-ifelse(IRF1_IRF2_DKO$adj.P.Val < 0.0005, 0.0005, 
                                IRF1_IRF2_DKO$adj.P.Val)

IRF1_IRF2_DKO$coef<-factor(IRF1_IRF2_DKO$coef,
                           levels = c("Ut_IRF1KO","Ut_IRF2KO","Ut_IRF1KO:IRF2KO"))  
ggplot(IRF1_IRF2_DKO,aes(x= coef,y = ensg, 
                         size = -log10(adj.P.Val),
                         color=logFC))+
  geom_point()+
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text = element_text(size = 12)) +
  ggtitle("Untreated")
ggsave(paste0(dir,"ut_IRF1_or_IRF2KO_TEST.png"),w=7,h=length(genes)*0.25+1)
#for other treatments
# Define the list of treatments
treatments <- c("IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
unique(limmaRes$genotype)
for (treatment1 in unique(limmaRes$treatment)){
  top_genes <- limmaRes %>%
    filter(treatment == treatment1, 
           genotype %in% c("IRF1KO","IRF2KO","IRF12KO","IRFDKO-SKOs"))%>%
    group_by(coef)%>%
    filter(adj.P.Val< 0.05) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = 20)%>%
    pull(ensg)%>%unique()
  
  # Filter limmaRes for the relevant genes and coefficients
  plot_data <- limmaRes %>%
    filter(ensg %in% top_genes,
           genotype %in% c("IRF1KO", "IRF2KO", "IRF12KO", "IRFDKO-SKOs","IRFDKO-SKOs"),
           treatment == treatment1)
  
  # Adjust the adj.P.Val for plotting
  plot_data$adj.P.Val <- ifelse(plot_data$adj.P.Val < 0.0005, 0.0005,
                                plot_data$adj.P.Val)
  
  plot_data$coef<-gsub(paste0(treatment1),"",plot_data$coef)
  plot_data$coef<-gsub("_|:","",plot_data$coef)
  plot_data$coef<-gsub("IRF1KOIRF2KO","IRF12KO",plot_data$coef)
  
  plot_data$coef<-factor(plot_data$coef,
                         levels = c("IRF1KO","IRF2KO","IRF12KO","IRFDKO-SKOs"))
  
  # Generate the plot
  ggplot(plot_data, aes(x = coef, y = ensg, 
                        size = -log10(adj.P.Val),
                        color = logFC)) +
    geom_point() +
    scale_color_gradient2(high = "red", low = "blue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle(paste0(treatment1))
  ggsave(paste0(dir,treatment1,"Top_20_genes.png"),h=length(top_genes)*0.25+1,w=7)
  
  
  dat.list <- list()
  for(gg in top_genes){
    dat.list[[gg]] <- metadata %>%
      mutate(E=scale(dataVoom$E[gg,])) %>%
      rownames_to_column("samples") %>%
      remove_rownames()
  }
  p.vals <- bind_rows(dat.list, .id="ensg") 
  
  plot<- p.vals%>%
    filter(treatment_time %in% c("Ut",treatment1))
  # Define the desired order of samples
  
  
  unique(plot$samples)
  # Convert 'sample' to a factor with the defined order
  plot$samples <- factor(plot$samples,
                         levels = c("Ut_Rosa_1",
                                    "Ut_Rosa_2",
                                    "Ut_Rosa_3",
                                    paste0(treatment1,"_Rosa_1"),
                                    paste0(treatment1,"_Rosa_2"),
                                    paste0(treatment1,"_Rosa_3"),
                                    "Ut_IRF1_1",
                                    "Ut_IRF1_2",
                                    "Ut_IRF1_3",
                                    paste0(treatment1,"_IRF1_1"),
                                    paste0(treatment1,"_IRF1_2"),
                                    paste0(treatment1,"_IRF1_3"),
                                    "Ut_IRF2_1",
                                    "Ut_IRF2_2",
                                    "Ut_IRF2_3",
                                    paste0(treatment1,"_IRF2_1"),
                                    paste0(treatment1,"_IRF2_2"),
                                    paste0(treatment1,"_IRF2_3"),
                                    "Ut_IRF12_1",
                                    "Ut_IRF12_2",
                                    "Ut_IRF12_3",
                                    paste0(treatment1,"_IRF12_1"),
                                    paste0(treatment1,"_IRF12_2"),
                                    paste0(treatment1,"_IRF12_3")))
  
  
  
  #mutate(stimulus = as.character(stimulus)) %>%
  ggplot(plot, aes(x=samples, y=ensg, fill=E)) + 
    geom_tile() +
    #facet_grid(. ~ stimulus, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text = element_text(size = 12)) 
  
  ggsave(paste0(dir,treatment1,"expression.png"),h=length(top_genes)*0.25+2,
         w=10)
}
}
##############################

dat.list <- list()
for(gg in genes){
  dat.list[[gg]] <- metadata %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
}
p.vals <- bind_rows(dat.list, .id="ensg") 

plot<- p.vals%>%
  filter(treatment_time %in% "Ut")
# Define the desired order of samples

# Convert 'sample' to a factor with the defined order
plot$samples <-factor(plot$samples,
                      levels = c("Ut_Rosa_1",
                                 "Ut_Rosa_2","Ut_Rosa_3",
                                 "Ut_IRF1_1","Ut_IRF1_2","Ut_IRF1_3",
                                 "Ut_IRF2_1","Ut_IRF2_2","Ut_IRF2_3",
                                 "Ut_IRF12_1","Ut_IRF12_2","Ut_IRF12_3"))


#mutate(stimulus = as.character(stimulus)) %>%
ggplot(plot, aes(x=samples, y=ensg, fill=E)) + 
  geom_tile() +
  #facet_grid(. ~ stimulus, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  theme(axis.text = element_text(size = 8)) 

ggsave("ut_IRF1_or_IRF2KO_TESTexpression.png")
########
#clustering
# Perform hierarchical clustering on 'E'
# Define the specific treatments and genotypes
treatments <- c("IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
genotypes <- c("IRF1KO", "IRF2KO", "IRF1KO:IRF2KO")

# Create a list to store the top genes for each combination of treatment and genotype
# Create a list to store the top genes for each combination of treatment and genotype




# Loop through each genotype and treatment


for (treatment1 in unique(limmaRes$treatment)){
  plot_genes <- limmaRes %>%
    filter(treatment == treatment1)%>%
    filter(genotype %in% c("IRF1KO","IRF2KO","IRF12KO","IRF2KO_vs_IRF1KO"))%>%
    group_by(coef)%>%
    filter(adj.P.Val< 0.05) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = 20)%>%
    ungroup%>%
    pull(ensg)%>%
    unique()
  
  plot_data <- limmaRes %>%
    filter(treatment == treatment1, genotype %in%  c("IRF1KO","IRF2KO","IRF12KO","IRF2KO_vs_IRF1KO"),
           ensg %in% plot_genes)
  
  # Plotting
  ggplot(plot_data, aes(x = reorder(ensg, -logFC), y = logFC, fill = color)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") + # Add black border
    facet_grid(rows = vars(genotype), scales = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(
      x = "Gene",
      y = "logFC",
      title = paste("Top 20 Genes by logFC for", treatment1, "Treatment"),
      fill = "Adj. P-Value"
    ) +
    scale_fill_manual(values = c("Significant" = "purple", "N.S" = "grey")) +
    guides(fill = guide_legend(title = "Adj. P-Value"))
  
  # Save each plot with a unique name based on treatment1
  ggsave(paste0(treatment1, ".png"), width = 25)
}


###############################################################
#
top_genes<-limmaRes %>%
  filter(genotype)
group_by(coef)%>%
  filter(adj.P.Val< 0.05) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 1000)%>%
  ungroup%>%
  pull(ensg)%>%unique()



# look at the coefficient names
colnames(coef(limmaFit))

# make sure we have the right names, otherwise we have to adapt the next line
stopifnot(all(colnames(coef(limmaFit)) == c("(Intercept)", "stimulusIFNa", "organSpleen", "stimulusIFNa:organSpleen")))

# now create a contrast matrix
contrast.mt <- cbind(IFNa_Spleen = c(0,1,0,1)) # we add the 2nd and 4th coefficient.
row.names(contrast.mt) <- colnames(coef(limmaFit))

# look at the matrix
contrast.mt

# Contrast fit similar to the original limma fit
limmaFit.contrast <- contrasts.fit(limmaFit,contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)

# Extract results for this contrast coefficient
limmaRes.contrast <- topTable(limmaFit.contrast, coef=colnames(contrast.mt),number = Inf) |>
  rownames_to_column("ensg") |>
  mutate(coef=colnames(contrast.mt))

# add them to the full table
limmaRes <- rbind(limmaRes.contrast, limmaRes) # add this coefficient to the result table
table(limmaRes$coef)
unique(limmaRes$coef)
###############################
#plots -----------------------
###############################
get_top_genes <- function(data, coefs_of_interest=unique(data$coef), top_n) {
  data %>%
    filter(coef %in% coefs_of_interest) %>%
    group_by(coef) %>%
    filter(adj.P.Val < 0.05) %>%
    filter(abs(logFC) > 1)%>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

coefs_of_interest = grep("IRF2KO_vs_IRF1KO",unique(limmaRes$coef),value = T)

IRF2_vs_IRF1 <- get_top_genes(limmaRes, coefs_of_interest, 200)
IRF2_vs_IRF1%>%
  arrange(abs(logFC))

IRF2_vs_IRF1_genes<-IRF2_vs_IRF1$ensg%>%unique()
data_plot<-dataVoom$E[IRF2_vs_IRF1_genes,
                      grep("Rosa",colnames(dataVoom$E),value = T)]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
colnames(data_plot)
#
####################################################
breaksList<-seq(-3,3,by =0.5)
col<- hcl.colors(11, "RdYlBu")
col1<-rev(col)# reverse the order of selected colours
#z-score normalization
data_subset_norm <- t(apply(data_plot, 1,cal_z_score))
colnames(data_subset_norm)
#####################################################
data_subset_norm<-data_subset_norm[,c("Ut_Rosa_1","Ut_Rosa_2","Ut_Rosa_3",
                                      "IFNb_4h_Rosa_1","IFNb_4h_Rosa_2","IFNb_4h_Rosa_3",
                                      "IFNb_24h_Rosa_1","IFNb_24h_Rosa_2","IFNb_24h_Rosa_3",
                                      "IFNg_4h_Rosa_1","IFNg_4h_Rosa_2","IFNg_4h_Rosa_3",
                                      "IFNg_24h_Rosa_1","IFNg_24h_Rosa_2","IFNg_24h_Rosa_3")]

# Identify duplicated rows (keep the first occurrence, remove the rest)
data_subset_norm <- data_subset_norm[!duplicated(data_subset_norm), ]

# Perform hierarchical clustering
row_clusters <- hclust(dist(data_subset_norm), method = "ward.D2")

# Cut the tree into 5 clusters
cluster_assignments <- cutree(row_clusters, k = 4)

# Get the order of rows as per the clustering
ordered_rows <- row_clusters$order

# Reorder clusters from top to bottom
ordered_cluster_assignments <- cluster_assignments[ordered_rows]

# Map the original cluster labels to ordered cluster labels from 1 to 5
unique_clusters <- unique(ordered_cluster_assignments)
cluster_order <- order(unique_clusters)
new_cluster_labels <- setNames(cluster_order, unique_clusters)
final_cluster_assignments <- new_cluster_labels[ordered_cluster_assignments]

# Create a data frame for row annotations with ordered cluster labels
annotation_row <- data.frame(Cluster = factor(final_cluster_assignments))
rownames(annotation_row) <- rownames(data_subset_norm)[ordered_rows]

# Generate the heatmap with annotations
png("test_clustering_IRF2_vs_IRF1_cluster_4_200.png")
pheatmap(data_subset_norm, 
         color = col1, 
         breaks = breaksList,
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         clustering_method = "ward.D2",
         show_colnames = TRUE, 
         show_rownames = FALSE,
         annotation_row = annotation_row,
         cutree_rows = 4)
dev.off()
genes_clusters <- data.frame(genes = rownames(data_subset_norm)[ordered_rows], 
                             cluster = final_cluster_assignments)
write.table(genes_clusters, paste0(dir,"genes_clusters_4_100.tsv"), row.names = FALSE)


dataVoom_IFNb<- dataVoom$E[,c(grep("Ut",colnames(dataVoom$E)),
                              grep("IFNb",colnames(dataVoom$E)))]
dataVoom_IFNb<-as.data.frame(dataVoom_IFNb)
dataVoom_IFNb <- dataVoom_IFNb %>% rownames_to_column(var = "genes")


dataVoom_long_IFNb <- dataVoom_IFNb %>% pivot_longer(cols =c(Ut_IRF1_1:IFNb_24h_Rosa_3), 
                                                     names_to = "sample", 
                                                     #rownames_to_column("genes"),
                                                     values_to = "cts",
                                                     values_transform = list(cts=as.numeric))
sample_info_IFNb<-metadata[c(grep("Ut",colnames(dataVoom$E)),grep("IFNb",
                                                                  colnames(dataVoom$E))),]
sample_info_IFNb["sample"]<-rownames(sample_info_IFNb)

dataVoom_long_IFNb <- full_join(dataVoom_long_IFNb, sample_info_IFNb, by = "sample")
# Create the genotype column based on the conditions
dataVoom_long_IFNb <- dataVoom_long_IFNb %>%
  mutate(genotype = case_when(
    IRF1 == "KO" & IRF2 == "KO" ~ "IRF12KO",
    IRF1 == "KO" ~ "IRF1KO",
    IRF2 == "KO" ~ "IRF2KO",
    TRUE ~ "WT"
  ))

head(dataVoom_long_IFNb)
#######
#candidate_genes
candidate_genes <- rownames(data_subset_norm)
data_wt_ko_long_IFNb<- dataVoom_long_IFNb %>% filter(genes %in% candidate_genes) %>%
  group_by(genes)%>%
  mutate(cts_zscore=(cts-mean(cts))/sd(cts))%>%
  group_by(genes,genotype,treatment_time)%>%
  #take the mean within one genotype_timepoint across the three replicate for a specific gene)
  summarise(mean_zscore_of_replicates = mean(cts_zscore),
            nrep=n())%>%ungroup()


dataVoom_cluster_IFNb <- data_wt_ko_long_IFNb %>% 
  inner_join(genes_clusters, by = "genes")

dataVoom_cluster_IFNb$genotype
dataVoom_cluster_IFNb$genotype <- factor(dataVoom_cluster_IFNb$genotype,
                                         levels = c("WT", "IRF1KO","IRF2KO","IRF12KO"))
#################


dataVoom_cluster_IFNb %>% 
  ggplot(aes(treatment_time, mean_zscore_of_replicates)) +
  geom_line(aes(group = genes), alpha = 0.1, colour ="gray") +
  geom_line(stat = "summary", fun = "median", 
            colour = "cadetblue4", size = 0.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(genotype), cols = vars(cluster))+
  theme(axis.text=element_text(size=8,angle = 90),
        axis.title=element_text(size=10,face="bold"))
ggsave(paste0(dir,"cluster_IFNb_4_200_line.png"))
######
#genotype X treatment
#####
#boxplot

# Assuming dataVoom_cluster_IFNb is already defined

dataVoom_cluster_IFNb %>%
  ggplot(aes(x = genotype, y = mean_zscore_of_replicates)) +
  geom_boxplot(aes(group = genotype), outlier.shape = NA, alpha = 0.6, fill = "cadetblue4") +
  #geom_jitter(colour = "black", width = 0.2, alpha = 0.4) +
  facet_grid(rows = vars(treatment_time), cols = vars(cluster)) +
  theme(axis.text = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10, face = "bold")) +
  labs(title = "Boxplot with Jitter for Each Cluster and Genotype",
       x = "Treatment Time",
       y = "Mean Z-score of Replicates")
ggsave(paste0(dir,"IFNb_cluster_4_200_boxplot.png"))


dataVoom_cluster_IFNb %>%
  ggplot(aes(x = genotype, y = mean_zscore_of_replicates)) +
  geom_violin(aes(group = genotype), alpha = 0.6, fill = "cadetblue4") +
  stat_summary(fun = median, geom = "point", size = 3, color = "black", shape = 95) +
  facet_grid(rows = vars(treatment_time), cols = vars(cluster)) +
  theme(axis.text = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10, face = "bold")) +
  labs(title = "Violin Plot with Median Line for Each Cluster and Genotype",
       x = "Genotype",
       y = "Mean Z-score of Replicates")
ggsave(paste0(dir,"IFNb_cluster_4_200_violin.png"))

######
# treatment X genotype
#####


dataVoom_cluster_IFNb %>%
  ggplot(aes(x = treatment_time, y = mean_zscore_of_replicates)) +
  geom_boxplot(aes(group = treatment_time), outlier.shape = NA, alpha = 0.6, fill = "cadetblue4") +
  #geom_jitter(colour = "black", width = 0.2, alpha = 0.4) +
  facet_grid(rows = vars(genotype), cols = vars(cluster)) +
  theme(axis.text = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10, face = "bold")) +
  labs(title = "Boxplot with Jitter for Each Cluster and Genotype",
       x = "Treatment Time",
       y = "Mean Z-score of Replicates")
ggsave(paste0(dir,"IFNb_cluster_4_200_boxplot_reverse.png"))


dataVoom_cluster_IFNb %>%
  ggplot(aes(x = treatment_time, y = mean_zscore_of_replicates)) +
  geom_violin(aes(group = treatment_time), alpha = 0.6, fill = "cadetblue4") +
  stat_summary(fun = median, geom = "point", size = 3, color = "black", shape = 95) +
  facet_grid(rows = vars(genotype), cols = vars(cluster)) +
  theme(axis.text = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10, face = "bold")) +
  labs(title = "Violin Plot with Median Line for Each Cluster and Genotype",
       x = "Genotype",
       y = "Mean Z-score of Replicates")
ggsave(paste0(dir,"IFNb_cluster_4_200_violin_reverse.png"))

#########
#
#
#
#
for (treatment1 in unique(limmaRes$treatment)){
  top_genes <- limmaRes %>%
    filter(treatment == treatment1, 
           genotype %in% c("IRF12KO"))%>%
    group_by(coef)%>%
    filter(adj.P.Val< 0.05) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = 20)%>%
    pull(ensg)%>%unique()
  
  # Filter limmaRes for the relevant genes and coefficients
  plot_data <- limmaRes %>%
    filter(ensg %in% top_genes,
           genotype %in% c("IRF1KO", "IRF2KO", "IRF12KO"),
           treatment == treatment1)
  
  # Adjust the adj.P.Val for plotting
  plot_data$adj.P.Val <- ifelse(plot_data$adj.P.Val < 0.0005, 0.0005,
                                plot_data$adj.P.Val)
  
  plot_data$coef<-gsub(paste0(treatment1),"",plot_data$coef)
  plot_data$coef<-gsub("_|:","",plot_data$coef)
  plot_data$coef<-gsub("IRF1KOIRF2KO","IRF12KO",plot_data$coef)
  
  plot_data$coef<-factor(plot_data$coef,
                         levels = c("IRF1KO","IRF2KO","IRF12KO"))
  
  # Generate the plot
  ggplot(plot_data, aes(x = coef, y = ensg, 
                        size = -log10(adj.P.Val),
                        color = logFC)) +
    geom_point() +
    scale_color_gradient2(high = "red", low = "blue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle(paste0(treatment1))
  ggsave(paste0(dir,treatment1,"Top_20_genesTEST.png"),h=length(top_genes)*0.25+1,w=7)
  
  
  dat.list <- list()
  for(gg in top_genes){
    dat.list[[gg]] <- metadata %>%
      mutate(E=scale(dataVoom$E[gg,])) %>%
      rownames_to_column("samples") %>%
      remove_rownames()
  }
  p.vals <- bind_rows(dat.list, .id="ensg") 
  
  plot<- p.vals%>%
    filter(treatment_time %in% c("Ut",treatment1))
  # Define the desired order of samples
  
  
  unique(plot$samples)
  # Convert 'sample' to a factor with the defined orde
  dat.list <- list()
  for(gg in top_genes){
    dat.list[[gg]] <- metadata %>%
      mutate(E=scale(dataVoom$E[gg,])) %>%
      rownames_to_column("samples") %>%
      remove_rownames()
  }
  p.vals <- bind_rows(dat.list, .id="ensg") 
  
  plot<- p.vals%>%
    filter(treatment_time %in% "Ut")
  # Define the desired order of samples
  
  # Convert 'sample' to a factor with the defined order
  if (treatment1 == "Ut") {
    plot$samples <- factor(plot$samples,
                           levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                      "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                      "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                      "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3"))
  } else {
    plot$samples <- factor(plot$samples,
                           levels = c("Ut_Rosa_1", "Ut_Rosa_2", "Ut_Rosa_3",
                                      paste0(treatment1, "_Rosa_1"), paste0(treatment1, "_Rosa_2"), paste0(treatment1, "_Rosa_3"),
                                      "Ut_IRF1_1", "Ut_IRF1_2", "Ut_IRF1_3",
                                      paste0(treatment1, "_IRF1_1"), paste0(treatment1, "_IRF1_2"), paste0(treatment1, "_IRF1_3"),
                                      "Ut_IRF2_1", "Ut_IRF2_2", "Ut_IRF2_3",
                                      paste0(treatment1, "_IRF2_1"), paste0(treatment1, "_IRF2_2"), paste0(treatment1, "_IRF2_3"),
                                      "Ut_IRF12_1", "Ut_IRF12_2", "Ut_IRF12_3",
                                      paste0(treatment1, "_IRF12_1"), paste0(treatment1, "_IRF12_2"), paste0(treatment1, "_IRF12_3")))
  }
  
  
  #mutate(stimulus = as.character(stimulus)) %>%
  ggplot(plot, aes(x=samples, y=ensg, fill=E)) + 
    geom_tile() +
    #facet_grid(. ~ stimulus, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text = element_text(size = 12)) 
  
  ggsave(paste0(dir,treatment1,"expression_TEST.png"),h=length(top_genes)*0.30+2,
         w=10)}

}
##########
#
# Example design matrix with invalid column names
design <- model.matrix(~ IRF1 * IRF2 * treatment_time-1, data = metadata)
# Check current levels of treatment_time
levels(metadata$treatment_time)

# If "Ut" is missing, add it to the levels
metadata$treatment_time <- factor(metadata$treatment_time, levels = c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h"))
design <- model.matrix(~ treatment_time * IRF1 * IRF2, data = metadata)
colnames(design)
####
#
#