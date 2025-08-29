###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(msigdbr)
################################################################################
inDir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
out<-"Ag_ScRNA_10_Pseudobulk_ex_in_NTC_functional_annotation/"
outdir<-dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_functional_annotation/")

################################################################################
#load limma results
dataVoom_NTC<-read_rds(inDir("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC<-read_rds(inDir("limma_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex<-read_rds(inDir("NTC_meta.rds"))
################################################################################
# Create a list with separate entries for each celltype-group pair
library(dplyr)
library(purrr)

# Create and save separate data frames for each celltype-group pair
limmaRes_NTC %>%
  group_by(celltype, group) %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 500) %>%
  group_split(celltype, group) %>%
  walk(function(df) {
    # Generate a unique file name based on celltype and group
    celltype <- unique(df$celltype)
    group <- unique(df$group)
    file_name <- outdir(paste0("data_", celltype, "_", group, ".csv"))
    
    # Save each dataframe as a CSV file
    write.csv(df, file_name, row.names = FALSE)
  })

  
# # # Download gene sets ------------------------------------------------------
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
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
#save(enr.terms, file=out("Genesets_Mouse.RData"))

################################################################################
#Enrichr------------------------------------------------------------------------
################################################################################
perform_enrichment_analysis <- function(limmaRes, direction) {
  enr.res.list <- list()
  
  for (ct in unique(limmaRes$celltype)) {
    # Extract genes of interest (GOI) for a given coefficient and direction
    goi <- limmaRes %>%
      filter(group == direction & celltype == ct) %>%
      pull(genes)
    
    # Perform enrichment analysis
    enr.res <- enrichr(goi, databases = c("KEGG_2019_Mouse",
                                          "MSigDB_Hallmark_2020",
                                          "WikiPathways_2019_Mouse",
                                          "GO_Biological_Process_2021",
                                          "TRRUST_Transcription_Factors_2019"))
    
    # Combine results into one long table
    enr.res <- bind_rows(enr.res, .id = "db")
    
    # Store results in the list
    enr.res.list[[ct]] <- enr.res
  }
  
  # Combine results for all cell types
  enr.res.all <- bind_rows(enr.res.list, .id = "celltype")
  
  return(enr.res.all)
}

# Perform enrichment analysis for up-regulated genes
limmaRes_NTC<-limmaRes_NTC[limmaRes_NTC$celltype != "MEP",]
enr.res.all_up <- perform_enrichment_analysis(limmaRes_NTC, "up")
write_rds(enr.res.all_up, ENRICHR("enr.res.all_NTC_up.rds"))

# Perform enrichment analysis for down-regulated genes
enr.res.all_down <- perform_enrichment_analysis(limmaRes_NTC, "down")
write_rds(enr.res.all_down, ENRICHR("enr.res.all_NTC_down.rds"))
head(enr.res.all_down)
################################################################################
#plot
################################################################################
enr.res.all_up <- read_rds(ENRICHR("enr.res.all_NTC_up.rds"))
enr.res.all_down <- read_rds(ENRICHR("enr.res.all_NTC_down.rds"))
dfs<-list(enr.res.all_up,enr.res.all_down)
names(dfs)<-c("up","down")
databases = c("KEGG_2019_Mouse" ,
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse", 
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")


pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db & df$celltype != "MEP", ]
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = celltype, 
                               y = Term, color = log2(Odds.Ratio),
                               size = pmin(10,neg.log10.Adjusted.P.value))) +
        
        geom_point() +
        scale_size_continuous(
          range = c(2, 6)  # Set the sizes in the legend
          )+
        scale_color_gradientn(
          colors = c("pink", "red"),
          #breaks = c(0, 5),  # Set custom breaks for the color scale
          #labels = c("0", "2", "5"),  # Set custom labels for the breaks
          #limits = c(0, 10)  # Set the limits for the color scale
        ) +
      
        ggtitle(paste0(db))+
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size=15),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12),
              title = element_text(size=18))+
        
      
      ggsave(ENRICHR(paste0(name, "_per.celltype_", db, ".pdf")), 
             w = 10, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})

##################
#without filter
##################
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db, ]
    df_filtered <- df_filtered %>% mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        pull(Term)
      enr<-df_filtered %>%filter(Term %in% plotting_enr)
      
      ggplot(enr, aes(x = celltype, 
                               y = Term, color = log2(Odds.Ratio),
                               size = pmin(10,neg.log10.Adjusted.P.value))) +
        
        geom_point() +
        scale_size_continuous(
          range = c(2, 6)  # Set the sizes in the legend
        )+
        scale_color_gradientn(
          colors = c("pink", "red"),
          #breaks = c(0, 5),  # Set custom breaks for the color scale
          #labels = c("0", "2", "5"),  # Set custom labels for the breaks
          #limits = c(0, 10)  # Set the limits for the color scale
        ) 
      #facet_wrap(vars(db))
      
      ggsave(ENRICHR(paste0(name, "_per.celltype_", db, "without_filter.pdf")), 
             w = 10, h = length(unique(enr$Term)) * 0.2 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})
################################################################################
# fgsea ------------------------------------------------------------------------
################################################################################
gsea.res <- data.table()  # Initialize the result table
for (ct in unique(limmaRes_NTC$celltype)) {
  #for (de.grp in unique(limmaRes_NTC[limmaRes_NTC$celltype == ct, ]$coef)) {
    for (dbx in names(enr.terms)) {
      # limmaRes_NTC_subset <- limmaRes_NTC[limmaRes_NTC$celltype == ct, ]
      # Correct subsetting with both conditions
      subset_limmaRes_NTC <- limmaRes_NTC[limmaRes_NTC$celltype == ct, ]
      # Now using `with` correctly on the subsetted data
      stats <- with(subset_limmaRes_NTC, setNames(logFC, nm = genes))
      
      # Check for missing values in stats
      if (any(is.na(stats))) {
        next  # Skip this iteration if there are missing values in stats
      }
      
      # Perform fgsea analysis
      fgsea_output <- fgsea(
        pathways = enr.terms[[dbx]],
        stats = stats
        #minSize = 15,   # Example additional arguments, adjust as necessary
        #maxSize = 500,  # Example additional arguments, adjust as necessary
        #nperm = 1000    # Example additional arguments, adjust as necessary
      )
      
      # Check if fgsea output is not empty
      if (length(fgsea_output) > 0) {
        gsea.res <- rbind(gsea.res, data.table(fgsea_output,
                                               #coef = de.grp,
                                               celltype = ct,
                                               db = dbx))
      }
    }
  }
gsea.res %>% write_rds(basedir("NTC_fgsea.rds"))
################################################################################
#plot fgsea---------------------------------------------------------------------
################################################################################
gsea.res <- read_rds(basedir("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge,
                                      function(vec) paste(vec[1:10], collapse = ","))

for(dbx in unique(gsea.res$db)){
  dat<-dirout(paste0(base,"FGSEA/",dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_",dbx,".tsv"))
}

dbx<-"MSigDB_Hallmark_2020"
# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  
  pDT <- gsea.res[db == dbx]
  ## Splitting the task to handle both ends of the NES spectrum-positive and negative
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  pDT <- pDT[pathway %in% pw.display]
  if (nrow(pDT) > 0) {
    # Proceed with your operation
    pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
    ggplot(pDT, aes(x=celltype, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
      
      scale_color_gradient2(low="blue", mid="white", high="red") +
      geom_point() +
      scale_size_continuous(range=c(0,5), limits = c(0,5)) +
      theme_bw(12) +
      xRot() +
      #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
      labs(x = "celltype")+
      
      theme(axis.text = element_text(size = 10)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
    #theme(strip.text.y=element_text(angle=0))
    dat<-dirout(paste0(base,"FGSEA/",dbx))
    ggsave(dat("GSEA_plot_",dbx,".png"), w=20,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
  } else {
    warning("Data table is empty. Skipping operation.")
  }
  
  #pDT <- hierarch.ordering(pDT, "pathway", "celltype", "NES", TRUE)
  
}

################################
