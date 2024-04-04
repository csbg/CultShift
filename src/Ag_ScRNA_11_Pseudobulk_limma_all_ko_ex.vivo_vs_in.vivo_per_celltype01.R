####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#14-03-24
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
#####################################################################
inDir<-dirout("/Ag_ScRNA_08_Pseudobulk")
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/"
basedir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/")
#####################################################################
#source DE function
#####################################################################
source("src/Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_perCT_function.R")
################################################
#load data and clean metadata
################################################
#metadata
meta<-read.csv(inDir("metadata.csv"),row.names=1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","Gran.")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

#only select the genotypes present in both tissue conditions
meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]

# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
#counts
counts <- read.csv(inDir("ex_counts.csv"), row.names = 1)
counts<-cbind(counts,read.csv(inDir("in_counts.csv"), row.names = 1))
counts<-counts[,rownames(meta)]
stopifnot(all(colnames(counts)==rownames(meta)))
stopifnot(all(unique(meta[meta$tissue=="ex.vivo",]$celltype)==unique(meta[meta$tissue=="in.vivo",]$celltype)))

################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
################################################
#perform DE independent of celltype
model_formula <- "~tissue*genotype"
result <- performDE(meta, counts,model_formula)
limmaRes<-bind_rows(result,.id="celltype")
head(limmaRes)
dataVoom<-result$dataVoom
limmaRes%>%write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT.rds"))

limmaRes[limmaRes$coef=="Ash1l" & limmaRes$group=="down",]
#for each coef
# Get unique cell types
unique(limmaRes$celltype)


# Iterate over unique coefficients
for(coefx in unique(limmaRes$coef)) {
  data <- limmaRes[limmaRes$coef == coefx, ]
  
  # Plotting for each unique cell type
  #for(celltype in unique_celltypes) {
    #cell_data <- data[data$celltype == celltype, ]
    
    # Plotting
    p <- ggplot(data, aes(x = logFC, y = pmin(-log10(adj.P.Val),5), col = group)) +
      geom_point() +
      scale_color_manual(values = c("down" = "#5782A7", "n.s" = "#B9B8B6", "up" = "#8A264A"), drop = FALSE) +
      ggtitle(paste0(gsub("ex.vivo:", "interaction:", coefx), "D.E")) +
      facet_wrap(vars(celltype), scales = "free_x")+
      #theme with white background
      theme_bw() +
      
      #eliminates background, gridlines, and chart border
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
      ) +
      
      #draws x and y axis line
      theme(axis.line = element_line(color = 'black'))
      
    # Facet wrapping by cell type
    
    # Save the plot
    out <- dirout(paste0(base,"/per_cell_type/",coefx))  # Adjust dirout function as per your setup
    ggsave(filename = out(paste0(coefx,"DE.pdf")), plot = p)
}
#######################################################################################################
#
for(ct in unique(limmaRes$celltype)) {
  data <- limmaRes[limmaRes$celltype == ct, ]
  
  # Plotting for each unique cell type
  #for(celltype in unique_celltypes) {
  #cell_data <- data[data$celltype == celltype, ]
  
  # Plotting
  p <- ggplot(data, aes(x = logFC, y = pmin(-log10(adj.P.Val),5), col = group)) +
    geom_point() +
    scale_color_manual(values = c("down" = "#5782A7", "n.s" = "#B9B8B6", "up" = "#8A264A"), drop = FALSE) +
    ggtitle(paste0(ct, "_D.E")) +
    facet_wrap(vars(coef), scales = "free_x")+
    #theme with white background
    theme_bw() +
    
    #eliminates background, gridlines, and chart border
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #panel.border = element_blank()
    ) +
    
    #draws x and y axis line
    theme(axis.line = element_line(color = 'black'))
  
  # Facet wrapping by cell type
  
  # Save the plot
  out <- dirout(paste0(base,"/per_cell_type/",ct))  # Adjust dirout function as per your setup
  
  ggsave(filename = out(paste0(ct,"DE.pdf")), plot = p,w=40,h=30)
}
#
counts <- table(limmaRes$group, limmaRes$KO)

# Convert counts to data frame
library(ggplot2)

# Assuming limmaRes is your data frame containing the "group" column and KO identifiers

# Filter out "n.s" group
library(ggplot2)

# Assuming limmaRes is your data frame containing the "group", "KO", and "celltype" columns

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]
head(limmaRes)
# Count occurrences of "up" and "down" genes for each KO within each cell type
counts <- with(limmaRes_filtered, table(celltype, group, coef))

# Convert counts to data frame
counts_df <- as.data.frame(counts)
names(counts_df) <- c("celltype", "group", "coef", "count")

# Plot bar plots for each cell type
plots <- lapply(unique(counts_df$celltype), function(celltype) {
  ggplot(subset(counts_df, celltype == celltype), aes(x = coef, y = count, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("up" = "#8A264A", "down" = "#5782A7")) +
    labs(title = paste("Number of 'up' and 'down' genes for each KO -", celltype),
         x = "KO", y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
})
out<-dirout(paste0(base,"counts"))

# Assuming limmaRes is your data frame containing the "group", "KO", and "celltype" columns

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]

# Remove duplicates to avoid counting duplicates multiple times
limmaRes_filtered_unique <- unique(limmaRes_filtered[, c("group", "KO", "celltype")])

# Count occurrences of "up" and "down" genes for each KO within each cell type
counts <- with(limmaRes_filtered_unique, table(celltype, group, KO))

# Convert counts to data frame
counts_df <- as.data.frame(counts)
names(counts_df) <- c("celltype", "group", "KO", "count")

# Create separate bar plots for each cell type
library(ggplot2)
library(gridExtra)

# Assuming limmaRes is your data frame containing the "group", "KO", and "celltype" columns

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]

# Remove duplicates to avoid counting duplicates multiple times
limmaRes_filtered_unique <- unique(limmaRes_filtered[, c("group", "KO", "celltype")])

# Count occurrences of "up" and "down" genes for each KO within each cell type
counts <- with(limmaRes_filtered_unique, table(celltype, group, KO))

# Convert counts to data frame
counts_df <- as.data.frame(counts)
names(counts_df) <- c("celltype", "group", "KO", "count")

# Create separate bar plots for each cell type
plots <- lapply(unique(counts_df$celltype), function(celltype) {
  df <- subset(counts_df, celltype == celltype)
  
  ggplot(df, aes(x = KO, y = count, fill = group)) +
    geom_bar(stat = "identity", position = "dodge", width = 1.5) +  # Adjust width as needed
    scale_fill_manual(values = c("up" = "#8A264A", "down" = "#5782A7")) +
    labs(title = paste("Number of 'up' and 'down' genes for each KO -", celltype),
         x = "KO", y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
})

# Arrange plots in a grid
grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol as needed
#########################################################

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]

# Count occurrences of "up" and "down" genes for each coef within each cell type
counts <- with(limmaRes_filtered, table(celltype, coef, group))
counts$type<-
# Convert counts to data frame
counts_df <- as.data.frame(counts)
names(counts_df) <- c("celltype", "coef", "group", "count")

# Create bar plot with facets for each cell type
ggplot(counts_df, aes(x = coef, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("up" = "#8A264A", "down" = "#5782A7")) +
  labs(title = "Number of 'up' and 'down' genes for each coef",
       x = "Coefficient", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(rows=vars(celltype))

library(ggplot2)

# Assuming limmaRes is your data frame containing the "group", "KO", "celltype", and "coef" columns

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]

# Count occurrences of "up" and "down" genes for each coef within each cell type
counts <- with(limmaRes_filtered, table(celltype, coef, group))

# Convert counts to data frame
counts_df <- as.data.frame(counts)
names(counts_df) <- c("celltype", "coef", "group", "count")

# Create bar plot with facets for each cell type
library(ggplot2)

# Assuming limmaRes is your data frame containing the "group", "KO", "celltype", and "coef" columns


# Create a new column "type" based on the condition
limmaRes$type <- ifelse(grepl("^ex\\.vivo:", limmaRes$coef), "interaction",
                        ifelse(limmaRes$coef == "ex.vivo", "ex.vivo", "KO"))


# Assuming limmaRes is your data frame containing the "group", "KO", "celltype", and "coef" columns

# Filter out "n.s" group
limmaRes_filtered <- limmaRes[limmaRes$group != "n.s", ]

# Count occurrences of "up" and "down" genes for each coef within each cell type and group
counts <- limmaRes_filtered %>%
  group_by(celltype, coef, group) %>%
  summarise(count = n())

# Create a new column "type" based on the condition
counts$type <- ifelse(startsWith(counts$coef, "ex.vivo:"), "interaction",
                      ifelse(counts$coef == "ex.vivo", "ex.vivo", "KO"))

# Create bar plot with facets for each cell type
ggplot(counts, aes(x = coef, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("up" = "#8A264A", "down" = "#5782A7")) +
  labs(title = "Number of 'up' and 'down' genes for each coef",
       x = "Coefficient", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cols=vars(type),rows = vars(celltype), scales = "free_x") +
  guides(fill = guide_legend(title = "Group"))

################################################
#****************#perform DE per celltype*******------
################################################

# Use purrr functions to iterate over unique KO and celltypes
ct<-unique(meta$celltype)[1]

map(unique(meta$celltype), ~ {
  ct <- .x
  meta_ct <- meta[meta$celltype==ct,] 
  
  if(nrow(meta_ct[meta_ct$tissue=="in.vivo",])>1 &
     nrow(meta_ct[meta_ct$tissue=="ex.vivo",])>1){
    meta_ct$genotype <- factor(meta_ct$genotype,
                               levels=c("NTC",unique(setdiff(meta_ct$genotype,"NTC"))))
    
    meta_ct$tissue <- factor(meta_ct$tissue, levels=c("in.vivo", "ex.vivo"))
    
    # Call the function for differential expression analysis
    limmaRes_ct <- performDE(meta_ct, counts,model_formula)
    
    # Add celltype column to the combined results
    limmaRes_ct$celltype <- ct
    
    # Combine with the existing limmaRes data frame
    limmaRes <<- rbind(limmaRes, limmaRes_ct)
    limmaRes$group <- ifelse(limmaRes$logFC >= 1 &
                               limmaRes$adj.P.Val <= 0.05, "up",
                             ifelse(limmaRes$logFC <= -1 &
                                      limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
  }else{
    print("null")
  }
  limmaRes
})

#############################################################
#ex.vivo condition
#############################################################
#volcano plot
coefx<-"ex.vivo"
data<-limmaRes[limmaRes$coef=="ex.vivo",]
ggplot(data,aes(x=logFC,y=-log10(adj.P.Val),col=group))+
  geom_point()+
  scale_color_manual(values = c("#5782A7", "#B9B8B6", "#8A264A"))+
  ggtitle(paste0("ex.vivo","_D.E"))
# facet_wrap(vars(data$celltype),scale="free")

out <- dirout(paste0(base, coefx))
ggsave(out("DE.pdf"))

#enrichment

source("src/Ag_enrichment.R")

perform_enrichment_analysis(data=data,
                            logFC_cutoff_up= 1,
                            logFC_cutoff_down= -3,
                            databases=c("KEGG_2019_Mouse", 
                                        "MSigDB_Hallmark_2020", 
                                        "WikiPathways_2019_Mouse", 
                                        "GO_Biological_Process_2021"),
                            output_file_prefix="all_CT")


####################################################




# Assuming your data frame is called 'data'

# Filter data based on the condition abs(logFC) > 1
gene_counts <- data %>%
  filter(abs(logFC) > 1)%>%
  group_by(genes) %>%
  summarize(unique_coefs = n_distinct(coef))

# Filter genes where the number of unique coefficients is greater than 1
overlap_genes <- gene_counts %>%
  filter(unique_coefs > 20) %>%
  pull(genes)

# Resulting list of genes that overlap across multiple coefficients with abs(logFC) > 1
overlap_genes
head(data)

pheatmap(data[genes %in% overlap_genes,])
##################################################################################
#interaction
out <- dirout(paste0(base,"INTERACTION"))
data<-limmaRes[limmaRes$coef %in% grep("ex.vivo:",limmaRes$coef,value=T),]

# Filter data based on the condition
filtered_data <- data %>%
  filter(adj.P.Val < 0.01) %>%
  group_by(coef) %>%
  summarize(genes = list(unique(ifelse(abs(logFC) > 2, ensg, NA)))) %>%
  unnest(genes)%>%na.omit()
head(filtered_data)
# Count occurrences of each gene
gene_counts <- unnested_data %>%
  count(genes)

# Filter genes that occur in multiple coefficients
genes_overlap <- gene_counts %>%
  filter(n > 1) %>%
  pull(genes)

# Resulting list of genes that overlap between coefficients


goi.all<-data %>%
  filter(adj.P.Val<0.01) %>%
  group_by(coef) %>%
  top_n(cutoff,abs(logFC)) %>%
  pull(ensg)%>%
  unique()
length(goi.all)

#


#up
cutoff<-10
top_pos_list<-list()
for (coefx in unique(data$coef)){
  top_pos<-data %>%
    filter(adj.P.Val<0.01 &coef==coefx) %>%
    slice_max(logFC,n=cutoff) %>%
    pull(ensg)%>%
    unique()
  top_pos_list[[coefx]]<-top_pos
}
goi.all<-unique(unlist(top_neg_list))
(p.coef <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) 
ggsave(out(paste0("top_pos_logFC_",cutoff,"_genes_across_KO.pdf")),
       w = 20, h = length(unique(goi.all)) * 0.2 + 3, limitsize = FALSE)
#top_down
top_neg_list<-list()
for (coefx in unique(data$coef)){
  top_neg<-data %>%
    filter(adj.P.Val<0.01 &coef==coefx) %>%
    slice_min(logFC,n=cutoff) %>%
    pull(ensg)%>%
    unique()
  top_neg_list[[coefx]]<-top_neg
}
goi.all<-unique(unlist(top_neg_list))
(p.coef <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) 
ggsave(out(paste0("top_neg_logFC_",cutoff,"_genes_across_KO.pdf")),
       w = 20, h = length(unique(goi.all)) * 0.2 + 3, limitsize = FALSE)

#####################################################################
#

selection<-"top_neg_logFC"
cutoff<-10
top_neg_list<-list()
for (coefx in unique(data$coef)){
  top_neg<-data %>%
    filter(adj.P.Val<0.01 &coef==coefx) %>%
    slice_min(logFC,n=cutoff) %>%
    pull(ensg)%>%
    unique()
  top_neg_list[[coefx]]<-top_neg
}
goi.all<-unique(unlist(top_neg_list))

(p.coef <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) 
ggsave(out(paste0("top_neg_logFC",cutoff,"_genes_across_KO.pdf")),
       w = 20, h = length(unique(goi.all)) * 0.2 + 3, limitsize = FALSE)

dat.list <- list()

for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
p.vals <- bind_rows(dat.list, .id="genes")
head(p.vals)
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ celltype, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
#interferon
selection<-"interferon_alpha_genes"
goi.all<-interferon_alpha_genes

# Find indices of empty character strings
goi.all <- Filter(nzchar, goi.all)
goi.all<-goi.all[goi.all %in% rownames(dataVoom)]
dat.list <- list()


(p.coef <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) 
ggsave(out(paste0("top_neg_logFC",cutoff,"_genes_across_KO.pdf")),
       w = 20, h = length(unique(goi.all)) * 0.2 + 3, limitsize = FALSE)

dat.list <- list()


for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
p.vals <- bind_rows(dat.list, .id="genes")

(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ celltype, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(rows = vars(celltype), space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 5)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample_cell_wise.pdf")),
       w = 40, h = 15)
unique(p.vals$sample1)
unique(p.vals$sample)
unique(p.vals$cell)
p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell))
dat<-ggplot(data=p.vals,aes(x=tissue, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ celltype, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dat1<-dat$data
ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
#####################################################################

load(PATHS$RESOURCES$Enrichr.mouse)
dbx<-unique(names(enr.terms))[1]


#pathway genes
inflammatory_response_genes<-enr.terms$MSigDB_Hallmark_2020$`Inflammatory Response`
interferon_alpha_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`
interferon_gamma_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
oxphos<-enr.terms$MSigDB_Hallmark_2020$`Oxidative Phosphorylation`
glycolysis<-enr.terms$MSigDB_Hallmark_2020$Glycolysis
random<-enr.terms$WikiPathways_2019_Mouse$`Dysregulated miRNA Targeting in Insulin/PI3K-AKT Signaling WP3855`
mtorc<-enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`
cholesterol<-enr.terms$WikiPathways_2019_Mouse$`Cholesterol Biosynthesis WP103`
adipogenesis<-enr.terms$WikiPathways_2019_Mouse$`Adipogenesis genes WP447`
fatty_acid<-enr.terms$MSigDB_Hallmark_2020$`Fatty Acid Metabolism`
Tnf<-enr.terms$MSigDB_Hallmark_2020$`TNF-alpha Signaling via NF-kB`
