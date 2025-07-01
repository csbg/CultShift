source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")
base<-"IRF2"
out <- dirout("IRF2")
countdata<-read.delim("/media/AGFORTELNY/PROJECTS/TfCf_AG/IRF2/raw_IRF2_gene.csv",sep = ",",header = T)
rownames(countdata) <- countdata$X
countdata$X<-NULL
nrow(countdata)


# Filtering
d0 <- DGEList(countdata)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
nrow(d)
######
#metadata
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
rownames(metadata) <- metadata$sample

# Ensure that the columns in countdata match the row names in metadata
countdata <- countdata[, rownames(metadata)]
stopifnot(all(colnames(countdata) == rownames(metadata)))
#

#edgeR norm
# Create the design matrix
design <- model.matrix(~treatment_time * IRF1 * IRF2, data = metadata)



png(out("voom.png"))
dataVoom<-voom(d, design, plot = TRUE)
dev.off()

############################################################
#fit model
colnames(design) <- make.names(colnames(design))
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)
colnames(limmaFit)
limmaRes <- list() # start an empty list
for(coefx in colnames(coef(limmaFit))){ # run a loop for each coef
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coef in a list.
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) %>%
    rownames_to_column("ensg")
}
limmaRes <- bind_rows(limmaRes, .id = "coef")
limmaRes$coef<-gsub("treatment_time","",limmaRes$coef)
limmaRes$coef<-gsub("genotype","",limmaRes$coef)
limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    coef %in% c("IRF1KO", "IRF2KO", "IRF1KO.IRF2KO") ~ paste0("Ut_", coef),
    TRUE ~ coef
  ))
unique(colnames(limmaFit$coefficients))
###############
#
#
# Assuming your design matrix is correctly set up with column names as shown before
# Sanitize column names to make them syntactically valid in R


# Check the new column names
##############################################

##########################################################
contrast.mt <- makeContrasts(
  Ut_doubleKO = IRF1KO + IRF2KO + IRF1KO.IRF2KO,
  levels = colnames(design)
)

# Fit the contrast
limmaFit.contrast <- contrasts.fit(limmaFit, contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)


limmaRes.contrast <- topTable(limmaFit.contrast, 
                              coef=colnames(contrast.mt),number = Inf) %>%
  rownames_to_column("ensg") %>%
  mutate(coef=colnames(contrast.mt))

# add them to the full table
limmaRes <- rbind(limmaRes.contrast, limmaRes) # add this coefficient to the result table

IRF2<-limmaRes%>% filter(coef=="Ut_IRF2KO")
##################################################
#
#
#

# Assuming `limmaRes` is already available and contains the necessary data
# Filter for relevant comparisons and reshape data
unique(limmaRes$coef)
logFC_data <- limmaRes %>%
  filter(coef %in% c(
    "Ut_IRF1KO","Ut_IRF2KO","Ut_doubleKO",
    grep("vs",unique(limmaRes$coef),value = T),
    grep("doubleKO",unique(limmaRes$coef),value = T)
    
  )) %>%
  select(ensg, coef, logFC) %>%
  pivot_wider(names_from = coef, values_from = logFC)
unique(logFC_data$)
# Remove rows with any NA values
logFC_data <- logFC_data %>%
  drop_na()

# Calculate pairwise correlations
cor_matrix <- cor(logFC_data[, -1], use = "pairwise.complete.obs")

# Create heatmap
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE)

###############################################
#fgsea
#############################################
source("src/00_init.R")
out <- dirout("EXT_02_EnrichR_Genesets/")
# # 
# # 
# # 
# # # Download gene sets ------------------------------------------------------
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
# 
# 
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
save(enr.terms, file=out("Genesets_Mouse.RData"))
##############################################
set.seed(123)  # For reproducibility
gsea.res <- data.table()
# Add a small amount of random noise to logFC to break ties
IRF2 <- IRF2 %>%
  mutate(logFC_no_tie = logFC + rnorm(n(), mean = 0, sd = 1e-5))


# Now using `with` correctly on the subsetted data
stats <- with(IRF2, setNames(logFC, nm = ensg))

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
for (dbx in names(enr.terms)) {
      # limmaRes_subset <- limmaRes[limmaRes$celltype == ct, ]
      # Correct subsetting with both conditions
      
      
      # Now using `with` correctly on the subsetted data
      stats <- with(IRF2, setNames(logFC, nm = ensg))
      
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
                                               
                                               
                                               db = dbx))
      }
    }

out<-dirout(paste0(base,"/FGSEA"))

saveRDS(gsea.res, file=out("FGSEA_IRF2KO.RDS"))
write_rds(enr.terms,basedir("mouse.enrichr.rds"))
# Plots -------------------------------------------------------------------
#if(!"gsea.res" %in% ls()) gsea.res <- readRDS(out("FGSEA.RDS"))
#load

gsea.res<-read_rds(out("FGSEA_IRF2KO.RDS"))
# cleanup / export results
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))

for(dbx in unique(gsea.res$db)){
  dat<-dirout(paste0(base,"FGSEA/",dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_",dbx,".tsv"))
}
head(gsea.res)
dbx<-"MSigDB_Hallmark_2020"
# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  
  pDT <- gsea.res[db == dbx]
  head()
  ## Splitting the task to handle both ends of the NES spectrum-positive and negative
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5)]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5)]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  pDT <- pDT[pathway %in% pw.display]
  #pDT <- hierarch.ordering(pDT, "pathway", "celltype", "NES", TRUE)
  pDT <- hierarch.ordering(pDT, "pathway", "NES",T)
  ggplot(pDT, aes(x=size, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    
    scale_color_gradient2(low="blue", mid="white", high="red") +
    #geom_point(data=pDT[padj < 0.05]) +
    scale_size_continuous(range=c(0,5), limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
    #labs(x = "")+
    
    theme(axis.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  #theme(strip.text.y=element_text(angle=0))
  dat<-dirout(paste0(base,"FGSEA/",dbx))
  ggsave(dat("GSEA_plot_",dbx,".pdf"), w=30,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
}


###############################################################################
#enrichr
###############################################################################
perform_enrichment_analysis <- function(limmaRes,
                                        databases,
                                        coefficients,
                                        logFC_threshold = 1,
                                        directory,
                                        file_prefix="") {
  enrichment_results_up <- list()
  enrichment_results_down <- list()
  
  # Iterate over each cell type and coefficient
  coefficients<-unique(limmaRes$coef)[unique(limmaRes$coef) %in% c("Ut_IRF2KO")]
  for (coefx in coefficients) {
    cat("Processing:", coefx, "\n")
    
    # Filter for upregulated and downregulated genes
    genes_up <- limmaRes %>%
      filter(coef == coefx & adj.P.Val <0.05 & logFC > logFC_threshold) %>%
      pull(ensg)
    
    genes_down <- limmaRes %>%
      filter(coef == coefx & adj.P.Val <0.05 & logFC < -logFC_threshold) %>%
      pull(ensg)
    
    if (length(genes_up) > 0) {
      enr_res_up <- enrichr(genes_up, databases = databases)
      if (!is.null(enr_res_up)) {
        # Filter out list elements with zero rows
        enr_res_up_filtered <- Filter(function(x) nrow(x) > 0, enr_res_up)
        if (length(enr_res_up_filtered) > 0) { # Check if there are any remaining list elements
          # Bind rows of the filtered list
          enr_res_up <- bind_rows(enr_res_up_filtered, .id = "db") %>%
            mutate(coef = coefx) # Add ct and coef as columns
          enrichment_results_up[[paste(coefx, "up", sep = "_")]] <- enr_res_up
        }
      }
    }
    
    if (length(genes_down) > 0) {
      enr_res_down <- enrichr(genes_down, databases = databases)
      if (!is.null(enr_res_down)) {
        # Filter out list elements with zero rows
        enr_res_down_filtered <- Filter(function(x) nrow(x) > 0, enr_res_down)
        if (length(enr_res_down_filtered) > 0) { # Check if there are any remaining list elements
          # Bind rows of the filtered list
          enr_res_down <- bind_rows(enr_res_down_filtered, .id = "db") %>%
            mutate(coef = coefx) # Add ct and coef as columns
          enrichment_results_down[[paste(coefx, "down", sep = "_")]] <- enr_res_down
        }
      }
    }
  }
  
  
  # Combine results
  
  down_enrichr <- bind_rows(enrichment_results_down, .id = "coef")
  up_enrichr <- bind_rows(enrichment_results_up, .id = "coef")
  
  #create directory for enrichr and save rds files
  
  #
  
  write_rds(down_enrichr, paste0(directory,"/",file_prefix,"down_logFC_",logFC_threshold,"_enrichr.rds"))
  write_rds(up_enrichr, paste0(directory,"/",file_prefix,"up_logFC_",logFC_threshold,"_enrichr.rds"))
  
}

perform_enrichment_analysis(
  limmaRes=limmaRes,
  databases="MSigDB_Hallmark_2020",
  coefficients=unique(limmaRes$coef),
  logFC_threshold = 1,
  directory = "/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/IRF2",
  file_prefix="")
enrichr_up<-read_rds("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/IRF2/up_logFC_1_enrichr.rds")
enrichr_down<-read_rds("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/IRF2/down_logFC_1_enrichr.rds")

################################################################################
enr.res<-enrichr_up
db="MSigDB_Hallmark_2020"
output_file = "Upregulated_enrichr"

  # Subset the enrichment results for the current database
  plotting_enr <- enr.res[enr.res$db == db,] %>%
    filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
    mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  # Check the number of overlapping genes
  plotting_enr <- plotting_enr %>%
    mutate(overlap_count = sapply(strsplit(Genes, ";"), length)) %>%
    filter(overlap_count > 1)
  
  # Create the plot
  ggplot(plotting_enr, aes(x = Overlap, 
                           y = Term, color = log2(Odds.Ratio),
                           size = neg.log10.Adjusted.P.value)) +
    labs(x = "IRF2KO",
         y = "Terms-Odds.Ratio > 5 & Adjusted.P.value < 0.05", 
         title = paste0(output_file))+
    geom_point() +
    
    # coord_flip()+
   # facet_wrap(vars(celltype))+
    theme_bw(12)+
    scale_size(range = c(1,3))+
    scale_color_gradient2(high="red", low="blue")+
    #scale_size_continuous(1,4)+
    ggtitle(paste0(output_file,"_",db))+
    theme(axis.title = element_text(size=7)) +
    theme(strip.text = element_text(size = 7))+
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 7),
          axis.text.y = element_text(size = 7))
  
  
  # Save the plot
  dir<-dirout(paste0(base,"/Enrichr/"))
  ggsave(filename = dir(paste0(output_file,"_", db, ".pdf")),
         #width = length(unique(coefficients))*0.35+4,
         #h=18,
         height = (length(unique(plotting_enr$Term)) * 0.155 + 3),
         limitsize = FALSE)
################################################################################
#
  enr.res<-enrichr_down
  db="MSigDB_Hallmark_2020"
  output_file = "Downregulated_enrichr"
  
  # Subset the enrichment results for the current database
  plotting_enr <- enr.res[enr.res$db == db,] %>%
    filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
    mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  # Check the number of overlapping genes
  plotting_enr <- plotting_enr %>%
    mutate(overlap_count = sapply(strsplit(Genes, ";"), length)) %>%
    filter(overlap_count > 1)
  
  # Create the plot
  ggplot(plotting_enr, aes(x = Overlap, 
                           y = Term, color = log2(Odds.Ratio),
                           size = neg.log10.Adjusted.P.value)) +
    labs(x = "IRF2KO",
         y = "Terms-Odds.Ratio > 5 & Adjusted.P.value < 0.05", 
         title = paste0(output_file))+
    geom_point() +
    
    # coord_flip()+
    # facet_wrap(vars(celltype))+
    theme_bw(12)+
    scale_size(range = c(1,3))+
    #scale_color_gradient2(high="red", low="blue")+
    #scale_size_continuous(1,4)+
    ggtitle(paste0(output_file,"_",db))+
    theme(axis.title = element_text(size=7)) +
    theme(strip.text = element_text(size = 7))+
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 7),
          axis.text.y = element_text(size = 7))
  
  
  # Save the plot
  dir<-dirout(paste0(base,"/Enrichr/"))
  ggsave(filename = dir(paste0(output_file,"_", db, ".pdf")),
         #width = length(unique(coefficients))*0.35+4,
         #h=18,
         height = (length(unique(plotting_enr$Term)) * 0.155 + 3),
         limitsize = FALSE)
  


























#function
extract_unique_genes <- function(enrichr_data,
                                 terms,
                                 db) {
  unique_genes <- enrichr_data %>%
    filter(db == db & Term %in% terms) %>%
    pull(Genes) %>%
    unique() %>%
    strsplit(";") %>%
    unlist() %>%
    toupper() %>%
    unique()
  return(unique_genes)
}
################################################################################

# Visualize results ---------------------------------------------------------
limmaRes <- res %>% filter(grepl("treatmentex_vivo",res$genotype))
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))


################################################################################
#function
perform_enrichment_analysis <- function(limmaRes,
                                        databases,
                                        coefficients,
                                        logFC_threshold = 1,
                                        directory,
                                        file_prefix="") {
  enrichment_results_up <- list()
  enrichment_results_down <- list() 
  
  # Iterate over each cell type and coefficient
  for (ct in unique(limmaRes$cell_type)) {
    cat("Processing:", ct, "\n")
    for (coefx in unique(coefficients)) {
      cat("Processing:", ct, coefx, "\n")
      
      # Filter for upregulated and downregulated genes
      genes_up <- limmaRes %>%
        filter(cell_type == ct & genotype == coefx & group == "up" & logFC > logFC_threshold) %>%
        pull(gene)
      
      genes_down <- limmaRes %>%
        filter(cell_type == ct & genotype == coefx & group == "down" & logFC < -logFC_threshold) %>%
        pull(gene)
      
      if (length(genes_up) > 0) {
        enr_res_up <- enrichr(genes_up, databases = databases)
        if (!is.null(enr_res_up)) {
          # Filter out list elements with zero rows
          enr_res_up_filtered <- Filter(function(x) nrow(x) > 0, enr_res_up)
          if (length(enr_res_up_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_up <- bind_rows(enr_res_up_filtered, .id = "db") %>%
              mutate(cell_type = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_up[[paste(ct, coefx, "up", sep = "_")]] <- enr_res_up
          }
        }
      }
      
      if (length(genes_down) > 0) {
        enr_res_down <- enrichr(genes_down, databases = databases)
        if (!is.null(enr_res_down)) {
          # Filter out list elements with zero rows
          enr_res_down_filtered <- Filter(function(x) nrow(x) > 0, enr_res_down)
          if (length(enr_res_down_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_down <- bind_rows(enr_res_down_filtered, .id = "db") %>%
              mutate(cell_type = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_down[[paste(ct, coefx, "down", sep = "_")]] <- enr_res_down
          }
        }
      }
    }
  }
  
  # Combine results 
  
  down_enrichr <- bind_rows(enrichment_results_down, .id = "ct_coef")
  up_enrichr <- bind_rows(enrichment_results_up, .id = "ct_coef")
  
  #create directory for enrichr and save rds files
  
  #
  
  write_rds(down_enrichr, paste0(directory,"/",file_prefix,"down_logFC_",logFC_threshold,"_enrichr.rds"))
  write_rds(up_enrichr, paste0(directory,"/",file_prefix,"up_logFC_",logFC_threshold,"_enrichr.rds"))
}
#function usage
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")
directory<-"/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr"

coefficients<-unique(limmaRes$genotype)
perform_enrichment_analysis(limmaRes, 
                            databases, coefficients,
                            directory = directory,logFC_threshold = 3)
################################################################################

################################################################################
IFN_limmRes <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Calculate percentage of downregulated genes (logFC < -1)
downregulated <- IFN_limmRes %>%
  filter(logFC < -1 & adj.P.Val < 0.05)
nrow(downregulated)
length(IFN_genes)
unique(IFN_genes)

percent_downregulated <- nrow(downregulated) / length(IFN_genes) * 100

# Create a bar plot for percentage of downregulated genes
percent_plot <- ggplot() +
  geom_bar(aes(x = "Downregulated", y = percent_downregulated), stat = "identity", fill = "blue") +
  geom_text(aes(x = "Downregulated", y = percent_downregulated + 2, label = paste0(round(percent_downregulated, 1), "%")), vjust = -0.5) +
  geom_bar(aes(x = "Not Downregulated", y = 100 - percent_downregulated), stat = "identity", fill = "gray") +
  geom_text(aes(x = "Not Downregulated", y = 100 - percent_downregulated + 2, label = paste0(round(100 - percent_downregulated, 1), "%")), vjust = -0.5) +
  coord_flip() +
  labs(x = "", y = "Percentage") +
  ggtitle("Percentage of Downregulated IFN Genes") +
  theme_minimal()

percent_plot
# Subset limmaRes for specific genes
limma_subset <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Filter limma results for IFN genes and treatmentex_vivo genotype
IFN_limmRes <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Count total IFN genes
total_IFN_genes <- length(IFN_genes)

# Count upregulated IFN genes (assuming logFC > 1 as upregulated)
upregulated_IFN <- IFN_limmRes %>%
  filter(logFC < -1) # Adjust threshold as needed

total_upregulated_IFN <- nrow(upregulated_IFN)

# Create a summary dataframe for plotting
summary_df <- data.frame(
  Category = c("Total IFN Genes", "Upregulated IFN Genes"),
  Count = c(total_IFN_genes, total_upregulated_IFN)
)
##########################
# Create the bar plot
 ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.3) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(
    title = "Summary of IFN Genes",
    x = "Category",
    y = "Count"
  ) +
  theme(legend.position = "none")
################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################
#for selected KOs
# plotting function for number of up/down genes
plot_genes <- function(data, 
                       count_threshold = 0,
                       width = 0.2,
                       title_suffix = "",
                       file_suffix = "") {
  
  # Filter out rows with count equal to 0 and based on count threshold
  filtered_data <- data %>% 
    filter(Count != 0) %>% 
    filter(Count >= count_threshold)
  
  # Plot
  p <- ggplot(filtered_data, 
              aes(x = genotype,
                  y = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) +
    geom_col(width = width) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    labs(x = "Interaction-KO",
         y = "Log10(Number of Genes)",
         title = paste("Upregulated and downregulated genes", title_suffix)) +
    theme_bw(12) +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_grid(cols =vars(cell_type)) +
    coord_flip() +
    theme(strip.text = element_text(size = 15))
  
  # Save plot
  ggsave(enrich(paste0("Number_of_genes", file_suffix, ".png")), plot = p,
         height = length(coefficients)*0.25+3)
}


#for all KOs
# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, genotype) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")

head(summary_df)
plot_genes(summary_df, count_threshold = 0, title_suffix = "", file_suffix = "")
#for Kos with above 10 genes differentially regulated between tissue type
plot_genes(summary_df, count_threshold = 10, title_suffix = " (Count >= 10)", 
           file_suffix = "_above10")

# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1


#filtered based on KOs with atleast 10 differentially regulated genes between tissue types
count_threshold = 10
coefficients <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(genotype)%>%
  unique()
perform_enrichment_analysis(limmaRes,
                            databases,
                            coefficients,
                            directory = directory,
                            file_prefix = "count_above10_")

################################################################################
#plot and save function for enriched results

plot_and_save<- function(enr.res, db, output_file,file_prefix = "",base) {
  # Subset the enrichment results for the current database
  plotting_enr <- enr.res[enr.res$db == db,] %>%
    filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
    mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  # Check the number of overlapping genes
  plotting_enr <- plotting_enr %>%
    mutate(overlap_count = sapply(strsplit(Genes, ";"), length)) %>%
    filter(overlap_count > 1)
  
  # Create the plot
  ggplot(plotting_enr, aes(x = gsub("ex.vivo:","",coef), 
                           y = Term, color = log2(Odds.Ratio),
                           size = neg.log10.Adjusted.P.value)) +
    labs(x = "KOs",
         y = "Terms-Odds.Ratio > 5 & Adjusted.P.value < 0.05", 
         title = paste0(output_file))+
    geom_point() +
    
    # coord_flip()+
    facet_wrap(vars(cell_type))+
    theme_bw(12)+
    scale_size(range = c(1,3))+
    #scale_size_continuous(1,4)+
    ggtitle(paste0(output_file,"_",db))+
    theme(axis.title = element_text(size=7)) +
    theme(strip.text = element_text(size = 7))+
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 7),
          axis.text.y = element_text(size = 7))
  
  
  # Save the plot
  enrichr<-dirout_jak(paste0("enrichr_plot/",db))
  
  ggsave(filename = enrichr(paste0(file_prefix,output_file,"_", db, ".pdf")),
         width = length(unique(coefficients))*0.35+4,
         #h=18,
         height = (length(unique(plotting_enr$Term)) * 0.155 + 3),
         limitsize = FALSE)
}



#set directory/filepath

enrich<-dirout_jak("Enrichr")
# Apply the function to each database using purrr::map()
down_enrichr<-read_rds(enrich("/down_logFC_3_enrichr.rds"))
up_enrichr<-read_rds(enrich("/up_logFC_3_enrichr.rds"))

map(databases, ~ plot_and_save(enr.res = up_enrichr, db = .x, output_file = "Upregulated_enrichr"))
map(databases, ~ plot_and_save(enr.res = down_enrichr, db = .x, output_file = "Downregulated_enrichr"))

################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################
#filteredKO
down_enrichr_filtered_KO<-read_rds(enrich("count_above10_down_logFC_1_enrichr.rds"))
up_enrichr_filtered_KO<-read_rds(enrich("count_above10_up_logFC_1_enrichr.rds"))
map(databases, ~ plot_and_save(enr.res = up_enrichr_filtered_KO, db = .x,
                               output_file = "Upregulated_enrichr",
                               file_prefix = "count_above10"))
map(databases, ~ plot_and_save(enr.res = down_enrichr_filtered_KO, db = .x,
                               output_file = "Downregulated_enrichr",
                               file_prefix = "count_above10"))

################################################################################
# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")
###############################################

# plotting function for number of up/down genes
plot_genes <- function(data, 
                       count_threshold = 0,
                       width = 0.2,
                       title_suffix = "",
                       file_suffix = "") {
  
  # Filter out rows with count equal to 0 and based on count threshold
  filtered_data <- data %>% 
    filter(Count != 0) %>% 
    filter(Count >= count_threshold)
  
  # Plot
  p <- ggplot(filtered_data, 
              aes(x = gsub("ex.vivo:", "", coef),
                  y = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) +
    geom_col(width = width) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    labs(x = "Interaction-KO",
         y = "Log10(Number of Genes)",
         title = paste("Upregulated and downregulated genes", title_suffix)) +
    theme_bw(12) +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_grid(cols =vars(cell_type)) +
    coord_flip() +
    theme(strip.text = element_text(size = 15))
  
  # Save plot
  ggsave(basedir(paste0("Number_of_genes", file_suffix, ".png")), plot = p,
         height = length(coefficients)*0.25+3)
}

#for all KOs
plot_genes(summary_df, count_threshold = 0, title_suffix = "", file_suffix = "")
#for Kos with above 10 genes differentially regulated between tissue type
plot_genes(summary_df, count_threshold = 10, title_suffix = "counts_above_10", 
           file_suffix = "_above10")
##########################
#do it by myself
# Filter out rows with count equal to 0 and based on count threshold
filtered_data <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)
head(summary_df)
# Plot
p <- ggplot(filtered_data, 
            aes(x = genotype,
                y = ifelse(Regulation == "Downregulated",
                           -log10(Count), log10(Count)), 
                fill = Regulation)) +
  geom_col(width = width) +
  scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
  labs(x = "Interaction-KO",
       y = "Log10(Number of Genes)",
       title = paste("Upregulated and downregulated genes", "counts_above_10")) +
  theme_bw(12) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cols =vars(cell_type)) +
  coord_flip() +
  theme(strip.text = element_text(size = 15))

# Save plot
basedir
ggsave(paste0(base,"Number_of_genes", "counts_above_10", ".pdf"), plot = p,
       )
################################################################################
#fgsea -------------------------------------------------------------------------
################################################################################
out <- dirout("EXT_02_EnrichR_Genesets/")
# # # Download gene sets ------------------------------------------------------
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
save(enr.terms, file=out("Genesets_Mouse.RData"))
################################################################################
################################################################################
#IFN-genes
################################################################################
up_enrichr<-read_rds(paste0(directory,"/","up_logFC_1_enrichr.rds"))
down_enrichr<-read_rds(paste0(directory,"/","down_logFC_1_enrichr.rds"))

#from enrichment analysis-

IFN_down<-extract_unique_genes(enrichr_data = down_enrichr, term=c("Interferon Alpha Response",
                                                                   "Interferon Gamma Response"),
                               db="MSigDB_Hallmark_2020")
##################################################
#all_downregulated_IFN genes across condtions
limmaRes%>%
  filter(toupper(gene) %in% IFN_down)%>%
  filter(adj.P.Val < 0.05)%>%
  filter(logFC < -1)%>%
  ggplot(aes(x = genotype, y = gene, 
             size = pmax(4,-log10(adj.P.Val)),
             color = pmin(10,logFC))) +
  geom_point() +
  scale_color_gradient2(high = "red", low = "blue") +
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(cell_type), scales = "free_x") +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
  ggtitle(paste("ex.vivo-in.vivo_JAK-STAT")) 

# Create directory if it doesn't exist

dir<-dirout_jak(paste0("/IFN_genes"))

# Save the plot
ggsave(dir(paste0("IFN_genes_downregulated","_summary.pdf")),
       width = length(IFN_down) * 0.05 ,
       height = length(IFN_down) * 0.1 + 2)
#limitsize = FALSE)
##################################################
# Define a function to summarize and plot IFN gene data
plot_Gene_set_summary <- function(gene_set, limmaRes, output_dir,gene_set_name=names(gene_set)) {
  
  # Filter and get unique genes from limmaRes based on the provided gene_set
  Gene_set <- intersect(gene_set, unique(limmaRes$gene))
  
  # Summarize data
  plot_data <- limmaRes %>%
    filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
    mutate(Regulation = case_when(
      logFC > 0 ~ "Upregulated",
      logFC < -1 ~ "Downregulated",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Regulation)) %>%
    group_by(genotype, cell_type, Regulation) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = genotype, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(. ~ cell_type, scales = "free_x", space = "free") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    theme_minimal() +
    labs(
      title = paste("Count of Upregulated and Downregulated Genes by Genotype and Cell Type"),
      subtitle = gene_set_name,
      x = "Genotype",
      y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")
  
  # Save the plot with a unique filename based on gene_set_name
  output_filename <- file.path(output_dir, paste0("Gene_set_summary_", gene_set_name, ".pdf"))
  ggsave(output_filename, plot = p, width = 10, height = 6)
  
  return(output_filename)
}

# Directory to save plots
output_dir <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr/"

# Iterate through each gene set name in pathway_names and create plots
plot_files <- lapply(pathway_genes, function(pathway_genes) {
  plot_Gene_set_summary(pathway_genes, limmaRes, output_dir)
})

# Print the list of saved plot files
print(plot_files)
# Example pathway names
pathway_names <- list(
  "Cholesterol Homeostasis",
  "mTORC1 Signaling",
  "Interferon Alpha Response",
  "Interferon Gamma Response",
  "Inflammatory Response",
  "p53 Pathway",
  "Myc Targets V1",
  "IL−6/JAK/STAT3 Signaling",
  "Oxidative Phosphorylation",
  "TGF−beta Signaling",
  "Unfolded Protein Response",
  "E2F Targets",
  "Protein Secretion"
)


# Initialize a list to store genes for each pathway



# Define function to summarize and plot gene set data
plot_Gene_set_summary <- function(gene_set_name, genes, limmaRes, output_dir) {
  
  # Extract genes from gene set
  Gene_set <- genes[[gene_set_name]]
  
  # Filter and get unique genes from limmaRes based on the provided gene_set
  Gene_set <- intersect(Gene_set, unique(limmaRes$gene))
  
  # Check if Gene_set is empty
  if (length(Gene_set) == 0) {
    message(paste("No genes found for", gene_set_name))
    return(NULL)
  }
  
  # Summarize data for upregulated and downregulated genes
  plot_data <- limmaRes %>%
    filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
    mutate(Regulation = case_when(
      logFC > 0 ~ "Upregulated",
      logFC < -1 ~ "Downregulated",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Regulation)) %>%
    group_by(genotype, cell_type, Regulation) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Check if plot_data is empty
  if (nrow(plot_data) == 0) {
    message(paste("No data found for", gene_set_name))
    return(NULL)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = genotype, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(. ~ cell_type, scales = "free_x", space = "free") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    theme_minimal() +
    labs(
      title = paste("Count of Upregulated and Downregulated Genes by Genotype and Cell Type"),
      subtitle = gene_set_name,
      x = "Genotype",
      y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")
  
  # Save the plot with a unique filename based on gene_set_name
  output_filename <- file.path(output_dir, paste0("Gene_set_summary_", gene_set_name, ".pdf"))
  ggsave(output_filename, plot = p, width = 10, height = 6)
  
  return(output_filename)
}

# Iterate through each gene set in pathway_genes and create plots
plot_files <- lapply(names(pathway_genes), function(gene_set_name) {
  plot_Gene_set_summary(gene_set_name, pathway_genes, limmaRes, output_dir)
})
length(ISG_core)
# Print the list of saved plot files
print(plot_files)

plot_data <- limmaRes %>%
  filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
  mutate(Regulation = case_when(
    logFC > 0 ~ "Upregulated",
    logFC < -1 ~ "Downregulated",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Regulation)) %>%
  group_by(genotype, cell_type, Regulation) %>%
  summarise(Count = n(), .groups = 'drop')


#####################
(tx <- SANN$tissue[1])
for(tx in unique(SANN$tissue)){
  out <- dirout(paste0(basedir, tx, "/"))
  monocle.file <- out("MonocleObject.RData")
  
  monocle.obj.list <- list()
  
  sx <- SANN[tissue == tx]$sample[1]
  for(sx in SANN[tissue == tx]$sample){
    fx <- inDir("SeuratObj_", sx, ".RData")
    if(!file.exists(fx)) stop(fx, " seurat object not found")
    
    print(paste("Reading ",sx))
    load(fx)
    
    # Add SoupX processing here
    soupChannel <- SoupChannel(rawCounts(seurat.obj), seurat.obj@assays$RNA@counts)
    soupChannel <- autoEstCont(soupChannel)
    cleaned <- adjustCounts(soupChannel)
    
    seurat.obj@assays$RNA@counts <- cleaned
    
    # Add more annotation
    for(x in c("tissue", "markers", "timepoint", "sample", "sample_broad")){
      seurat.obj@meta.data[[x]] <- SANN[sample == sx][[x]]
    }
    
    # Extract information for Monocle
    mat.use <- seurat.obj@assays$RNA@counts
    stopifnot(!any(duplicated(row.names(mat.use))))
    if(!"GFP" %in% row.names(mat.use)){
      x <- matrix(0, nrow = 2, ncol = ncol(mat.use))
      row.names(x) <- c("GFP", "BFP")
      mat.use <- rbind(mat.use, x)
    }
    monocle.obj.list[[sx]] <- new_cell_data_set(expression_data = mat.use, cell_metadata = seurat.obj@meta.data)
    
    # Store CITE-seq data
    if(sx == "DM_CITEseq-2_NA_NM_1"){
      citeseq.MT <- additional.info.x
      save(citeseq.MT, file=out.base("CITESEQ_Antibodies.RData"))
    }
  }
  
  # Make sure all objects have the same number of rows
  stopifnot(length(unique(sapply(monocle.obj.list, nrow))) == 1)
  
  # Combine datasets
  monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = FALSE)
  
  # Process dataset
  monocle.obj <-
    preprocess_cds(monocle.obj, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
  set.seed(42)
  monocle.obj <- align_cds(monocle.obj, 
                           alignment_group = "sample_broad", 
                           residual_model_formula_str = "~Phase",
                           verbose = TRUE
  )
  monocle.obj <- reduce_dimension(monocle.obj,
                                  reduction_method = "UMAP",
                                  preprocess_method = "Aligned",
                                  verbose = TRUE)
  
  # Clustering
  set.seed(12121)
  monocle.obj = cluster_cells(monocle.obj)
  
  # Store full dataset
  save(monocle.obj, file=monocle.file)
}
#############################
IRF2_down<-limmaRes%>%filter(coef=="Ut_IRF2KO")%>%
  filter(logFC< -1,adj.P.Val< 0.05)
dat.list<-list()
for(gg in unique(synergy_genes)) {
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
  
  ggsave(dir(paste0(gene,"_synergy.pdf")))
}
colnames(dat.list)
#For each KO where there were significantly up genes, make the plot
for (gg in synergy_genes){
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
pdf(out("heatmap_synergy_genes_genes.pdf"),h=16)
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Set to FALSE to avoid clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Adjust color scheme
         main = "Expression Heatmap",
         fontsize = 8)
dev.off()

