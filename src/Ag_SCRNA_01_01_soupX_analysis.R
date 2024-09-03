source("src/00_init.R")
library(SoupX)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
require("sceasy")
library(enrichR)
###############################################
Indir <- dirout("SCRNA_01_01_Seurat/")
out <- dirout("Ag_SCRNA_01_01_soupX_analysis/")

# Convert to mouse --------------------------------------------------------
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)

hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})
interferon_gamma_genes <- enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`
interferon_beta_genes <- enr.terms$MSigDB_Hallmark_2020$`Interferon Beta Response`
cholesterol_genes <- enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
# Combine both gene sets (if you want to consider overlaps between both responses)
interferon_genes <- union(interferon_gamma_genes, interferon_beta_genes)
hemoglobin_genes_mouse <- c("Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-y", "Hbb-bh1", "Hbb-bh2", "Hbz")

# Here counts- count in soupx pool, est count estimated in soup/ totalcount

# Write the final combined dataframe to a file

Annotation <- fread(Indir("SampleAnnotation.tsv"), fill=TRUE, sep="\t", header = TRUE)
soupx <- fread(Indir("combined_soupx_results_with_genes_and_est.tsv"),fill=TRUE, sep="\t", header = TRUE)
colnames(soupx)<-gsub("Sample","sample",colnames(soupx))
soupx_annotation <- inner_join(Annotation,soupx, by ="sample")
# Exclude rows where the Gene column is empty or NA
soupx_annotation <- soupx_annotation[!(soupx_annotation$Gene == "" | is.na(soupx_annotation$Gene)), ]
length(interferon_genes)


soupx_annotation <- soupx_annotation[,c("sample","tissue","timepoint","est","counts","Gene","Category","Total_Est")]

# Check unique values in 'tissue' to understand data spread
unique(soupx_annotation$tissue)
soupx_annotation <- soupx_annotation[tissue %in% c("in vivo","ex vivo")]
########################################################################
# Identify the top 10 genes by 'counts'
top_genes <- soupx_annotation %>%
  arrange(desc(counts)) %>%  # Sort by counts in descending order
  head(20) %>%               # Select top 10 rows
  pull(Gene)                 # Extract the Gene column

# Scatter plot of est counts for each gene, colored by tissue type
ggplot(soupx_annotation, aes(x = counts, y = est, shape = Category, color = tissue)) +
  geom_point(size = 3) +  # Make the points a bit larger
  #facet_wrap(~ tissue) +    # Create a separate plot for each gene
  theme_minimal() +       # Use a minimal theme
  theme(axis.text.x = NULL) + 
  xlim(0,20000)+
  ylim(0,0.005)+# Rotate x-axis labels for better readability
  labs(
    x = "Counts",
    y = "Estimated Soupx Counts",
    title = "Estimated Counts by Tissue Type for Each Gene",
    color = "Tissue"
  )+
  geom_text(
    data = soupx_annotation %>% filter(Gene %in% top_genes),  # Filter only top 10 genes
    aes(label = Gene),        # Add labels to points for top 10 genes
    vjust = -1,               # Position label above the point
    size = 3,                 # Size of the text
    check_overlap = TRUE      # Avoid overlapping text labels
  )
ggsave(out("Soupx_estimate.pdf"))
sc$soupProfile["Idi1",]
sum(seurat.obj@assays$RNA@counts["Bst2",])
test<-soupx_annotation %>% filter(counts > 100)%>%
  filter(est > 0.0001)%>%
  pull(Gene)%>%
  unique()
genes_to_exclude <- soupx_annotation %>% filter(counts > 100)%>%
  filter(est > 0.0005)%>%
  pull(Gene)%>%
  unique()
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
