####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#limma with all 
#data
#
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
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype/"
basedir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype/")
#####################################################################
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
#source DE function
#####################################################################
source("src/Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_function.R")
################################################
#out <- dirout(paste0(base,"/per_cell_type/"))
model_formula <- "~tissue*genotype"
limmaRes <- performDE(meta, counts, model_formula)
limmaRes<-limmaRes%>%
  filter(coef != "(Intercept)")
limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "genotype", "")) %>%
  mutate(coef = str_replace(coef, "tissue", ""))
unique(limmaRes$coef)
#only interaction
limmaRes<-limmaRes[limmaRes$coef %in% grep("ex.vivo:",limmaRes$coef,value=T),]%>%na.omit()
#dataVoom<-result$dataVoom
limmaRes%>%write_rds(basedir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
###########################################################################


#limmaRes[limmaRes$coef=="Ash1l" & limmaRes$group=="down",]
#for each coef
#############
#Volcano plots
################

# # Iterate over unique coefficients
# for(coefx in unique(limmaRes$coef)) {
#   data <- limmaRes[limmaRes$coef == coefx, ]
#   
#   # Plotting for each unique cell type
#   #for(celltype in unique_celltypes) {
#     #cell_data <- data[data$celltype == celltype, ]
#     
#     # Plotting
#     p <- ggplot(data, aes(x = logFC, y = pmin(-log10(adj.P.Val),5), col = group)) +
#       geom_point() +
#       scale_color_manual(values = c("down" = "#5782A7", "n.s" = "#B9B8B6", "up" = "#8A264A"), drop = FALSE) +
#       ggtitle(paste0(gsub("ex.vivo:", "interaction:", coefx), "D.E")) +
#       facet_wrap(vars(celltype), scales = "free_x")+
#       #theme with white background
#       theme_bw() +
#       
#       #eliminates background, gridlines, and chart border
#       theme(
#         plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank()
#       ) +
#       
#       #draws x and y axis line
#       theme(axis.line = element_line(color = 'black'))
#       
#     # Facet wrapping by cell type
#     
#     # Save the plot
#     out <- dirout(paste0(base,"/per_cell_type/",coefx))  # Adjust dirout function as per your setup
#     ggsave(filename = out(paste0(coefx,"DE.pdf")), plot = p)
# }

#############################################################
# #ex.vivo condition
# #############################################################
# #volcano plot
# 
# ggplot(data,aes(x=logFC,y=pmin(-log10(adj.P.Val),5),col=group))+
#   geom_point()+
#   facet_wrap(vars(celltype))+
#   scale_color_manual(values = c("#5782A7", "#B9B8B6", "#8A264A"))+
#   ggtitle(paste0("ex.vivo","_D.E"))
# # facet_wrap(vars(data$celltype),scale="free")
# 
# out <- dirout(paste0(base, coefx))
# ggsave(out("DE.pdf"))
# 
# #enrichment

# source("src/Ag_enrichment.R")
# 
# perform_enrichment_analysis(data=data,
#                             logFC_cutoff_up= 1,
#                             logFC_cutoff_down= -3,
#                             databases=c("KEGG_2019_Mouse", 
#                                         "MSigDB_Hallmark_2020", 
#                                         "WikiPathways_2019_Mouse", 
#                                         "GO_Biological_Process_2021"),
#                             output_file_prefix="all_CT")
# # 


####################################################
