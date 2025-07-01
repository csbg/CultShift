##################################################
#***Ag_ScRNA_13_invivo_izzo_limma********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#19-03-24
#Aarathy
###########
########################################################################
source("src/00_init.R")
inDir<-dirout("/Ag_ScRNA_08_Pseudobulk_limma")

base<-"Ag_ScRNA_09_Pseudobulk_limma_and_enrich_ex_izzo_NTC/"
basedir<-dirout("Ag_ScRNA_09_Pseudobulk_limma_and_enrich_ex_izzo_NTC/")
library(edgeR)
library(limma)
########################################################################
#load data and clean metadata
########################################################################
#metadata
meta<-read.delim(inDir("metadata.tsv"),row.names=1)
meta<-meta%>%filter(meta$tissue=="ex.vivo")
meta<-meta[,c("sample","celltype","tissue")]
NTC_meta<-meta[grep("NTC",rownames(meta),value = T),]
#rownames(NTC_meta)<-gsub("NA","",rownames(NTC_meta))

#izzo
meta1<-read.delim(inDir("izzo_metadata.tsv"),row.names=1)
meta1$sample<-gsub(".*(WT[0-9]+).*", "\\1", rownames(meta1))
meta1$tissue<-"izzo"

meta<-rbind(NTC_meta,meta1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T.cell","IMP2","IMP1","Neu","Gran.")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]


# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)\\-]", ".", rownames(meta))
#rownames(meta)<-gsub("..",".",rownames(meta))
# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))
rownames(meta) <- gsub("E/B", "E.B", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta[] <- lapply(meta, gsub, pattern = "E/B", replacement = "E.B.")
meta[] <- lapply(meta, gsub, pattern = "E.B.", replacement = "Eo.Ba")
unique(meta[meta$tissue=="ex.vivo",]$celltype)
unique(meta[meta$tissue=="izzo",]$celltype)

meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta<-meta[!grepl("T.cell",rownames(meta)),]

meta[meta$celltype=="Eo",]$celltype<-gsub("Eo","Eo.Ba",meta[meta$celltype=="Eo",]$celltype)
meta[meta$celltype=="Ba",]$celltype<-gsub("Ba","Eo.Ba",meta[meta$celltype=="Ba",]$celltype)


##############

#counts
counts_izzo<-read.csv(inDir("izzo_counts.csv"), row.names = 1)
counts_izzo<-counts[,rownames(meta[meta$tissue=="izzo,"])]%>%na.omit()
counts <-read.delim(inDir("combined_in_ex_counts.tsv"), row.names = 1)
counts_NTC<-counts[rownames(counts_izzo),rownames(meta[meta$tissue=="ex.vivo",])]
counts<-cbind(counts_NTC,counts_izzo)
counts<-counts[,rownames(meta)]%>%na.omit()

stopifnot(all(colnames(counts)==rownames(meta)))


unique(meta$tissue)
# Assuming you have two count matrices: count_matrix1 and count_matrix2
# Perform data preprocessing and normalization for each dataset separately

izzo_counts <- DGEList(counts = counts_izzo)

d0 <- calcNormFactors(izzo_counts)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

dge1 <- calcNormFactors(d)

dge2 <- DGEList(counts = count_matrix2)
dge2 <- calcNormFactors(dge2)

# Merge normalized count matrices
merged_counts <- merge(dge1$counts, dge2$counts, by = "row.names", all = TRUE)

# Perform additional normalization or batch effect correction if necessary

# Proceed with limma

#########################################
#Normalization in combined dataset
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group <- interaction(meta$celltype, meta$tissue)
unique(group)

mm <- model.matrix(~0 + group)
rownames(mm) <- rownames(meta)

# Normalization
dataVoom <- voom(d, mm)
dataVoom%>% write_rds(basedir("dataVoom_ex.vivo_vs_izzo.rds"))
#


#Expression Upregulated
##############
#NOTE:CHECK IF CAPITALIZED
gene_list<-IFN_genes
#Theses are the KO with significant upregulated interaction terms
list_ko <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

dat.list <- list()
IFN_genes<-limmaRes%>%filter(toupper(ensg) %in% IFN_genes)%>%pull(ensg)%>%unique
# Iterate over each unique gene
for(gg in unique(IFN_genes)) {
  # Subset the metadata and E values for the current gene
  gene_data <- meta %>%
    mutate(E = scale(dataVoom$E["Ifit3",])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  avg_gene_data <- gene_data %>%
    group_by(tissue, celltype) %>%
    summarise(avg_E = mean(E))%>%
    mutate(gene=gg)
  
  # Store the average gene data in the list
  dat.list[[gg]] <- avg_gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
head(gene_data)
ggplot(gene_data,aes(x=tissue, y=E)) + 
  geom_boxplot() +
  geom_jitter()+
  facet_grid(~ celltype, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))
ggsave(NTC_Expr(paste0(terms_of_interest,"_NTC_averagedsample.png")),w=10,h=8)




#Here we have sign.upregulated genes for each KO in all celltypes(not only the one with sig.up) of that KO.
limma_KO_up_IFN_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
limma_KO_up_IFN_genes%>%write_rds(basedir("IFN_genes_in_each_ko_upregulated_for all_celltypes"))
head(limma_KO_up_IFN_genes)
######################################################
# IFN genes-------------------------------------------
######################################################
IFN_UP<-dirout(paste0(base,"IFN_inter/UP/"))
dirout_load(IFN_UP)
# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (comp in unique(limma_KO_up_IFN_genes$comparison)){
  data<-limma_KO_up_IFN_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
#MODIFY FROM HERE
################################################
#factors and levels
################################################
#meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("izzo", "ex.vivo"))
################################################
#perform DE independent of celltype
model_formula <- "~tissue"
source("src/Ag_ScRNA_13_invivo_izzo_limma_function.R")
limmaRes <- performDE(meta, counts,model_formula="~tissue",plot_title = "ex.vivo-izzo")
limmaRes%>%
  write_rds(out("limma_ex_izzo_NTC_all_CT.rds"))
#per celltype
celltypes_to_exclude <- c("IMP1","IMP2","Neu","MEP","GMP")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

source("src/Ag_ScRNA_09_pseudobulk_per_celltype_limma_function.R")
#
######
performDE <- function(meta, counts, tissue_type1, tissue_type2) {
  # Filtering
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  
  group <- interaction(meta$celltype, meta$tissue)
  unique(group)
  
  mm <- model.matrix(~0 + group)
  rownames(mm) <- rownames(meta)
  
  # Normalization
  dataVoom <- voom(d, mm)
  limmaFit <- lmFit(dataVoom, mm)
  limmaFit <- eBayes(limmaFit)
  
  # Generate comparisons based on tissue types
  
  comparisons <- lapply(unique(meta$celltype), function(celltype) {
    comp1 <- paste0("group", celltype, ".", tissue_type1)
    comp2 <- paste0("group", celltype, ".", tissue_type2)
    list(comp1, comp2, celltype)
  })
  
  # Initialize list to store top table results
  top_table <- list()
  
  # Perform differential expression analysis for each comparison
  for (comp_list in comparisons) {
    group1 <- comp_list[[1]]
    group2 <- comp_list[[2]]
    id <- comp_list[[3]]
    
    # Create contrast matrix
    contrast <- paste0(group1, "-", group2)
    
    # Debugging: Print contrast matrix
    print(contrast)
    
    # Ensure that contrast inputs are accepted
    if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
      stop("Invalid column names for comparison:", group1, " and ", group2)
    }
    
    # Fit model and perform differential expression analysis
    tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
    tmp <- eBayes(tmp)
    top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
  }
  
  # Bind results and add metadata
  top_table_res <- bind_rows(top_table, .id = 'celltype') %>%
    mutate(genes = gsub("\\...\\d+", "", rownames(.)))
  
  # Create a column 'group' to group the genes as up, down, or ns
  top_table_res$group <- ifelse(top_table_res$logFC >= 1 & 
                                  top_table_res$adj.P.Val <= 0.05, "up", 
                                ifelse(top_table_res$logFC <= -1 & 
                                         top_table_res$adj.P.Val <= 0.05, "down", "n.s"))
  
  # Plotting and saving
  gg <- ggplot(data = top_table_res, aes(x = logFC, y = pmin(10, -log10(adj.P.Val)), col = group)) +
    geom_point() +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    ggtitle(paste0(tissue_type1, "_", tissue_type2, "_NTC_D.E")) +
    facet_wrap(vars(celltype), scales = "free_x") +
    theme_bw()
  out<-dirout(paste0(base,"/perCT"))
  ggsave(out(paste0(tissue_type1, "_", tissue_type2, "_NTC_D.E.pdf")), plot = gg,w=20,h=15)
  
  return(top_table_res)
}

######
########
ex_izzo_NTC_per_ct<-performDE(meta, counts, 
                            tissue_type1="ex.vivo", 
                            tissue_type2="izzo")
out <- dirout(paste0(base, "per_cell_type/"))
ex_izzo_NTC_per_ct %>%
  write_rds(out("limma_ex_izzo_NTCperCT.rds"))

#################
##################################
#enrichment function using enrichr
#################################
out<-dirout("Ag_ScRNA_13_exvivo_izzo_limma_and_enrich/")

source("src/Ag_enrichment.R")
perform_enrichment_analysis(data=limmaRes,
                            logFC_cutoff_up= 1,
                            logFC_cutoff_down= -1,
                            databases = c("KEGG_2019_Mouse",
                                          "MSigDB_Hallmark_2020",
                                          "WikiPathways_2019_Mouse",
                                          "GO_Biological_Process_2021"),
                            output_file_prefix="Enrichment_exvivo-izzo_all_celltypes")
