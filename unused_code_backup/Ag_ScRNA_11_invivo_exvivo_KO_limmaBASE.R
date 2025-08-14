##################################################
#***Ag_ScRNA_08_Limma_on_Pseudobulk_NTC.R********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#06-12-23
#Aarathy
#############
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(dplyr)
library(tidyverse)
library(enrichR)
library(purrr)
library(data.table)
(load("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis//EXT_02_EnrichR_Genesets/Genesets_Mouse.RData"))
PATHS
inDir<-dirout("/ScRNA_08_Pseudobulk_limma")
enrichment_per_ko<-"Ag_ScRNA_09_enrichment_per_ko/"

out<-dirout("Ag_ScRNA_09_Pseudobulk_limma_all_ko/")

########################
#load data
########################
#meta
meta<-read.delim(inDir("metadata.tsv"),row.names=1)

meta<-meta[meta$celltype!="B-cell" & meta$celltype!="CLP" & meta$celltype!="Ery" & meta$celltype!="EBMP" & meta$celltype!="unclear",]
meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]

#change rownames to match column names of expression data
rownames(meta)<-gsub("Eo/Ba","Eo.Ba",rownames(meta))
meta$cell<-gsub("Eo/Ba","Eo.Ba",meta$cell)
meta$celltype<-gsub("Eo/Ba","Eo.Ba",meta$celltype)
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
rownames(meta)<-gsub(" ",".",rownames(meta))
rownames(meta)<-gsub("\\(",".",rownames(meta))
rownames(meta)<-gsub("\\)",".",rownames(meta))
rownames(meta)<-gsub("-",".",rownames(meta))
##############
#factors and levels

meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
#counts
counts <- read.delim(inDir("combined_in_ex_counts.tsv"), row.names = 1)

counts<-counts[,rownames(meta)]
stopifnot(all(colnames(counts)==rownames(meta)))

##################
# Create a function for differential expression analysis
performDE <- function(meta, counts) {
  # edgeR based normalization
  d0 <- counts[, rownames(meta)]
  d0 <- DGEList(d0)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  
  # model matrix
  
  des <- model.matrix(~tissue*genotype, data=meta)
  
  #A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
  #A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
  #The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs
  # Define contrasts for each knockout
  # voom
  dataVoom <- voom(d, des, plot = T)
  limmaFit <- lmFit(dataVoom, des)
  limmaFit <- eBayes(limmaFit)
  
  
  # limma Results
  #limmaRes <- list() # start an empty list
  unique(colnames(limmaFit$coef))
  head(limmaFit)
  colnames(coef(limmaFit))
  limmaRes <- list() # start an empty list
  
  limmaRes <- map_dfr(unique(colnames(limmaFit$coef)), ~{
    coefx <- .x
    topTable(limmaFit, coef = coefx, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = coefx)
  }) %>%
    filter(coef != "(Intercept)")
  limmaRes <- limmaRes %>%
    mutate(coef = str_replace(coef, "genotype", "")) %>%
    #mutate(coef = str_replace(coef, "tissueex.vivo:", "interaction_")) %>%
    mutate(coef = str_replace(coef, "tissue", ""))
    #mutate(coef = str_replace(coef, "exvivo:", "interaction_")) 
} 
#NTC ex.vivo---------------------------------
###############################################################
#celltype-specific#
###############################################################
# List to store results for each KO per celltype
results_list <- list()
# Create an empty data frame to store combined results

limmaRes<-data.frame()
# Use purrr functions to iterate over unique KO and celltypes
ct<-unique(meta$celltype)[1]

limmaRes <- performDE(meta, counts)

map(unique(meta$celltype), ~ {
  ct <- .x
  meta_ct <- meta[meta$celltype==ct,] 
  
  if(nrow(meta_ct[meta_ct$tissue=="in.vivo",])>1 &
     nrow(meta_ct[meta_ct$tissue=="ex.vivo",])>1){
  meta_ct$genotype <-  factor(meta_ct$genotype, levels=c("NTC", unique(setdiff(meta_ct$genotype,"NTC"))))
 
  meta_ct$tissue <- factor(meta_ct$tissue, levels=c("in.vivo", "ex.vivo"))
    
      # Call the function for differential expression analysis
  limmaRes_ct <- performDE(meta_ct, counts)
  head(limmaRes)
      # Add KO and celltype columns to the combined results
  
  limmaRes_ct$celltype <- ct
      
      # Combine with the existing combined_results data frame
  combined_results <<- rbind(combined_results, limmaRes_ct)
  combined_results$group
      
      # Add results to the sublist
  results_list[[ct]] <- limmaRes
  }else{
  print("null")
  }
})
  
#############################################################
#enrichment
#############################################################
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                                limmaRes$adj.P.Val <= 0.05, "up", 
                                 ifelse(limmaRes$logFC <= -1 & 
                                          limmaRes$adj.P.Val <= 0.05, "down", "n.s"))

data=limmaRes[limmaRes$coef=="ex.vivo",]

ggplot(data,aes(x=logFC,y=-log10(adj.P.Val),col=group))+
  geom_point()+
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3"))+
  ggtitle("NTC_differential_expression-exvivo_vs_invivo")
 # facet_wrap(vars(data$celltype),scale="free")
ggsave(out("NTC_differential_expression_test_onlyNTC.pdf"))
##################################################################
#enrichment
##################################################################
#all celltypes
###############################
data=limmaRes[limmaRes$coef=="ex.vivo",]
perform_enrichment_analysis <- function(goi, databases, output_file) {
  enr.res <- enrichr(goi, databases = databases)
  enr.res <- bind_rows(enr.res, .id = "db")
  
  map(databases, ~{
    db <- .x
    plotting_enr <- enr.res[enr.res$db==db,] %>% 
      filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
      mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
    
    ggplot(plotting_enr, aes(x = db, y = Term, size = log2(Odds.Ratio),
                             color = neg.log10.Adjusted.P.value)) +
      geom_point() +
      #scale_color_gradient2(low = "blue", mid = "white", high ="red", 
                           #midpoint = median(plotting_enr$neg.log10.Adjusted.P.value), space = "rgb", guide = "colourbar") +
      #facet_grid(rows = vars(substr(db, 1, 10)), scales = "free_y") +  
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+theme_bw(12)
    
    ggsave(out(paste0(output_file, "_", db, ".pdf")), w=20,h=length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
    
  })
}


# Perform enrichment analysis for downregulated genes
goi_down <- data %>% filter(group == "down") %>% filter(logFC < -2) %>%
  pull(ensg) %>% unique()

perform_enrichment_analysis(goi_down, 
                            databases = c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                            output_file = "Downregulated_exvivo_allcelltype.pdf")

# Perform enrichment analysis for upregulated genes
goi_up <- data %>% filter(group == "up") %>% pull(ensg) %>% unique()

perform_enrichment_analysis(goi_up, 
                            databases = c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                            output_file = "upregulated_exvivo_allcelltype.pdf")

###############################
#celltype-specific
# Enrichment analysis
#########################
# Function to perform enrichment analysis for a given celltype
enr.res.up.list <- list()
enr.res.down.list <- list()

# Unique celltypes in combined_results
unique_celltypes <- unique(combined_results$celltype)

# Define a function to process each celltype
process_celltype <- function(ct) {
  df <- combined_results %>%
    filter(celltype == ct, coef == "ex.vivo")
  
  if (nrow(df) == 0) {
    warning(paste("No data for celltype:", ct))
    return(list(up = data.frame(), down = data.frame()))
  }
  
  goi_up <- df %>%
    filter(group == "up" & logFC > 1) %>%
    pull(ensg) %>%
    unique()
  
  goi_down <- df %>%
    filter(group == "down" & logFC < -1) %>%
    pull(ensg) %>%
    unique()
  
  enr.res_up <- enrichr(goi_up, databases = c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")) %>%
    bind_rows(.id = "db") %>%
    mutate(celltype = ct)
  
  enr.res_down <- enrichr(goi_down, databases = c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")) %>%
    bind_rows(.id = "db") %>%
    mutate(celltype = ct)
  
  return(list(up = enr.res_up, down = enr.res_down))
}

# Loop through unique celltypes and store results in lists
enr.res.list <- map(unique_celltypes, process_celltype)

# Combine the lists into data frames
up <- bind_rows(enr.res.list %>% map(pluck, "up"), .id = "ct")
down <- bind_rows(enr.res.list %>% map(pluck, "down"), .id = "ct")

databases<-c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")
dfs<-list(up,down)
names(dfs)<-c("up","down")
#plot
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db, ]
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = celltype, y = Term, size = log2(Odds.Ratio),
                               color = neg.log10.Adjusted.P.value)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
        theme_bw(12) +
        facet_wrap(vars(db))
      
      ggsave(out(paste0(name, "_per.celltype_", db, ".pdf")), 
             w = 20, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})

################################################################
#KO interaction-------------------------------------------------
################################################################
#celltype-specific
# Unique celltypes and KOs in combined_results
unique_celltypes <- unique(combined_results$celltype)
t<-grep("ex.vivo:",unique(combined_results$coef),value = T)
unique_kos <- unique(gsub("ex.vivo:", "", t))

# Extract interaction term
df <- combined_results %>%
  filter(coef %in% grep("ex.vivo:",coef,value = T))%>%
  mutate(coef=gsub("ex.vivo:","",coef))%>%
  na.omit()
head(df)
# Initialize lists to store results
enr.res.up.list <- list()
enr.res.down.list <- list()


# Loop through unique celltypes and KOs
for (ct in unique_celltypes) {
  for (KO in unique_kos) {
    df1 <- df %>%
      filter(coef == KO, celltype == ct)
    
    if (nrow(df1) == 0) {
      warning(paste("No data for KO:", KO, "and celltype:", ct))
      enr.res_up <- data.frame(db = character(), celltype = character(), KO = character())
      enr.res_down <- data.frame(db = character(), celltype = character(), KO = character())
    } else {
      goi_up <- df1 %>%
        filter(group == "up" & logFC > 1) %>%
        pull(ensg) %>%
        unique()
      
      goi_down <- df1 %>%
        filter(group == "down" & logFC < -1) %>%
        pull(ensg) %>%
        unique()
      
      enr.res_up <- if (length(goi_up) > 0) {
        tryCatch(
          enrichr(goi_up, databases = c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")) %>%
            bind_rows(.id = "db") %>%
            mutate(celltype = ct, KO = KO),
          error = function(e) {
            warning(paste("Error in enrichment for up-regulated genes. KO:", KO, "and celltype:", ct, "\n", e))
            data.frame(db = character(), celltype = character(), KO = character())
          }
        )
      } else {
        data.frame(db = character(), celltype = character(), KO = character())
      }
      
      enr.res_down <- if (length(goi_down) > 0) {
        tryCatch(
          enrichr(goi_down, databases = c("KEGG_2019_Mouse","MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")) %>%
            bind_rows(.id = "db") %>%
            mutate(celltype = ct, KO = KO),
          error = function(e) {
            warning(paste("Error in enrichment for down-regulated genes. KO:", KO, "and celltype:", ct, "\n", e))
            data.frame(db = character(), celltype = character(), KO = character())
          }
        )
      } else {
        data.frame(db = character(), celltype = character(), KO = character())
      }
    }
    
    enr.res_up_list_name <- paste(ct, KO, "up", sep = "_")
    enr.res_down_list_name <- paste(ct, KO, "down", sep = "_")
    
    enr.res.up.list[[enr.res_up_list_name]] <- enr.res_up
    enr.res.down.list[[enr.res_down_list_name]] <- enr.res_down
  }
}

# Combine the lists into data frames
# Initialize empty data frames
up <- data.frame()
down <- data.frame()

# Loop through lists and bind non-empty data frames
for (name in names(enr.res.up.list)) {
  df <- enr.res.up.list[[name]]
  if (nrow(df) > 0) {
    df$ct_KO <- rep(name, nrow(df))
    up <- rbind(up, df)
  }
}

for (name in names(enr.res.down.list)) {
  df <- enr.res.down.list[[name]]
  if (nrow(df) > 0) {
    df$ct_KO <- rep(name, nrow(df))
    down <- rbind(down, df)
  }
}
######################################################
#plot KO----------------------------------------------
######################################################


databases<-c("KEGG_2019_Mouse" ,"MSigDB_Hallmark_2020","WikiPathways_2019_Mouse", "GO_Biological_Process_2021")
dfs<-list(up,down)
names(dfs)<-c("up","down")
#plot
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db, ]
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = celltype, y = Term, size = log2(Odds.Ratio),
                               color = neg.log10.Adjusted.P.value)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
        theme_bw(12) +
        facet_grid(cols=vars(KO))
      
      ggsave(out(paste0(name, "_per.celltype_ko", db, ".pdf")), 
             w = 20, h = length(unique(plotting_enr$Term)) * 0.5 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})


#####################################################################################################################
#specigic genes and gene sets
###################################################################


goi.enr <- enr.res.up %>%
  filter(Adjusted.P.value < 0.01 & Odds.Ratio > 10) %>%
  pull(Genes) %>%
  strsplit(";") %>%
  unlist()%>%
  unique()
head(enriched)
enriched<-combined_results %>%
  #mutate(gene = gmap[ensg,]$external_gene_name) %>%
  filter(toupper(genes) %in% goi.enr)


#genes present in the enrichment and top genes logfc>2 adjp<0.01
set1<-enriched %>% 
  filter(adj.P.Val<0.01 & abs(logFC)>2)%>%pull(genes)
set1<-enriched[enriched$genes%in%set1,]

ggplot(set1, aes(x = id, y = genes, size = -log10(adj.P.Val), color = logFC)) +
  geom_point() +
  scale_color_gradient2(high = "red", low = "blue") +
  theme(text = element_text(size = 10)) 
ggsave(out("NTC_top_DEGs_from enriched_pathways.pdf"))
###################################################
goi.all<-combined_results %>%
  filter(adj.P.Val<0.01) %>%
  group_by(id) %>%
  top_n(10,logFC) %>%
  pull(genes)

(p.coef <- combined_results %>%
    filter(genes %in% goi.all) %>%
    ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw())
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
head(dataVoom$E)

(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
ggsave(out(""))
#interferon alpha genes

combined_results %>%
  filter(genes %in% interferon_alpha_genes) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("Interferon_alpha_reponse_genes_40.pdf"),w=5,h=7)

head(goi.all)
interferon_alpha_genes[34]
dat.list <- list()
for(gg in interferon_alpha_genes[10:40]){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
head(dataVoom$E)
dataVoom$E[,]
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))
#interferon gamma genes
combined_results %>%
  filter(genes %in% interferon_gamma_genes[1:50] ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("Interferon_gamma_reponse_genes_50.pdf"),w=5,h=7)
#inflammatory
combined_results %>%
  filter(genes %in% inflammatory_response_genes ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("Inflammatory_50.pdf"),w=5,h=7)
#cholesterol_homeostasis
combined_results %>%
  filter(genes %in% cholesterol_homestasis ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("cholesterol_homeostasis.pdf"),w=5,h=7)
#oxphos
combined_results %>%
  filter(genes %in% oxphos ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("Oxphos.pdf"),w=5,h=7)
#glycolisis
combined_results %>%
  filter(genes %in% glycolysis ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("glycolysis.pdf"),w=5,h=7)

#mtorc
combined_results %>%
  filter(genes %in% mtorc[1:150] ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("mtorc.pdf"),w=5,h=7)

#cholesterol
combined_results %>%
  filter(genes %in% cholesterol_homestasis ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("cholesterol_homeostasis.pdf"),w=5,h=7)

#adipogenesis
combined_results %>%
  filter(genes %in% adipogenesis ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("adipogenesis.pdf"),w=5,h=7)
#random
combined_results %>%
  filter(genes %in% random ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()
#fatty_acid
combined_results %>%
  filter(genes %in% fatty_acid ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("fatty_acid.pdf"),w=5,h=7)
#TNF
combined_results %>%
  filter(genes %in% Tnf ) %>%
  ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(out("Tnf.pdf"),w=5,h=7)
#######################
######################
#

dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#############################################
goi_log2_2<-set1<-enriched %>% 
  filter(adj.P.Val<0.01 & abs(logFC)>2)%>%
  pull(genes)

dat.list <- list()
for(gg in goi_log2_2){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=gsub("_*","",sample), y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
############################################
t<-combined_results %>% filter(adj.P.Val<0.01)%>%
  top_n(1,abs(logFC))%>%
  pull(genes)

head(t)
test2<-meta%>%
  mutate(exp=dataVoom$E["Gbp2",])
ggplot(test2,aes(x=tissue,y=exp))+
  geom_boxplot()+geom_jitter()+facet_grid(.~celltype)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(out("example_genes_Gbp2.pdf"))
########################################
#fgsea--------------------------------------------
###############################################
(load(PATHS$RESOURCES$Enrichr.mouse))
dbx<-unique(names(enr.terms))[1]
head(unique(names(enr.terms)))
for(ct in unique(combined_results$id)){
  for(dbx in names(enr.terms)){
    gsea.res <- rbind(gsea.res, data.table(fgsea(
      pathways=enr.terms[[dbx]], 
      stats=with(combined_results[id == ct], setNames(logFC, nm=genes))), 
      grp=ct,
      db=dbx))
    
  }
}

for(dbx in unique(gsea.res$db)){
  write.tsv(gsea.res.export[db == dbx], out("GSEA_significant_",dbx,".tsv"))
}


# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  pDT <- gsea.res[db == dbx]
  
  pw.display <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=10), by=c("grp")]$pathway)
  pDT <- pDT[pathway %in% pw.display]
  pDT <- hierarch.ordering(pDT, "pathway", "grp", "NES", TRUE)
  
  ggplot(pDT, aes(x=grp, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05], shape=1, color="black") +
    scale_size_continuous(range=c(0,5), limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    facet_grid(. ~ db, space="free", scales="free") +
    theme(strip.text.y=element_text(angle=0))
  ggsave(out("GSEA_plot_",dbx,".pdf"), w=20,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
}




saveRDS(gsea.res, file=out("FGSEA.RDS"))

# Plots -------------------------------------------------------------------
if(!"gsea.res" %in% ls()) gsea.res <- readRDS(out("FGSEA.RDS"))
gsea.res<-data.table(gsea.res)
# cleanup / export results
gsea.res[is.na(NES), NES := 0]

gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
head(gsea.res.export)
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, 
                                      function(vec) paste(vec[1:10], collapse = ","))

write.tsv(gsea.res.export, out("GSEA_significant_all",".tsv"))


for(dbx in unique(gsea.res.export$db)){
  write.tsv(gsea.res.export[db == dbx], out("GSEA_significant_all_",dbx,".tsv"))
}

all_gsea<-gsea.res.export
all_gsea<-all_gsea[db=="Cancer_Cell_Line_Encyclopedia"]
all_gsea<-all_gsea[,c("pathway","ES","grp")]

all_gsea<-all_gsea[!duplicated(all_gsea),]



all_gsea <- all_gsea%>%
  pivot_wider(names_from = grp, values_from = ES)

all_gsea<-as.data.frame(all_gsea)
rownames(all_gsea)<-all_gsea$pathway
all_gsea<-all_gsea[,-1]
head(all_gsea)
corMT <- as.matrix(cor(all_gsea,use = "pairwise.complete.obs"))
corMT[is.na(corMT)] <- 0
head(corMT)
library(ComplexHeatmap)


png("Cancer_Cell_Line_Encyclopedia.png")
Heatmap(corMT)
dev.off()

################################################################
#############################################
#genes_from_pathways
############################################
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
################################################
#-------Correlation of logfc between celltypes
#################################################################

split_data <- split(combined_results, combined_results$KO)

# Split the whole dataframe based on the 'coef' column


# Determine the unique order of genes across all dataframes
# Assuming you have a list of dataframes in 'splitdata'

# Define the order for sorting the "gene" column
all_genes <- unique(combined_results$ensg)
common_gene_order <- sort(all_genes)
head(common_gene_order)

# Define a function to sort genes within each coef across all dataframes
sort_genes_across_categories <- function(df_list) {
  # Iterate through each dataframe in the list
  lapply(df_list, function(df) {
    # Iterate through unique categories
    unique_categories <- unique(df$coef)
    
    # Sort genes within each coef
    df_sorted <- lapply(unique_categories, function(coef) {
      coef_df <- subset(df, coef == coef)
      coef_df$gene <- factor(coef_df$ensg, levels = common_gene_order, ordered = TRUE)
      coef_df <- coef_df[order(coef_df$ensg), ]
      return(coef_df)
    })
    
    # Combine the sorted dataframes for each coef
    do.call(rbind, df_sorted)
  })
}

# Apply the sorting function to the list of dataframes
split_data <- sort_genes_across_categories(split_data)
split_data<-lapply(split_data, function(df){
  df<-df[,c("ensg","celltype","logFC","KO")]
  df<-pivot_wider(df,id_cols=ensg,
                  names_from=celltype,
                  values_from=logFC)
  return(df)
})
head(split_data)
ordered_data<-lapply(ordered_data, function(df){
  df<-df[,c("id","genes","logFC")]
  df<-pivot_wider(df,id_cols=genes,
                  names_from=id,
                  values_from=logFC)
  return(df)
}

logFC_df<-ordered_data %>%
  reduce(full_join, by = "genes") %>%
  column_to_rownames(var = "genes")

corMT<-cor(logFC_df)
pdf(out("logFC_correlation_between_celltypes.pdf"))
pheatmap(corMT)
dev.off()
goi_up <- combined_results %>% filter(group=="up" & id==ct)
head(enr.res.up)
