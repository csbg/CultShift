source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)

# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
#out <- dir.create(out)
#out <- dirout(out)

# Settings (from JAK-STAT paper) ------------------------------------------
NORMALIZE <- TRUE
LMM.REML <- FALSE
LMM.VoomWeights <- FALSE
BPPARAM=BiocParallel::SnowParam(10, "SOCK", progressbar=TRUE)
CUTOFFS <- list(
  DE.ATAC = list(
    model="Lme4",
    AveExpr = 0,
    adj.P.Val = 0.05
  ),
  DE.RNA = list(
    model="Dream",
    AveExpr = 0,
    adj.P.Val = 0.05
  ),
  GSEA.REPEATS=10000
)

# Load data ---------------------------------------------------------------

# download and extract data
if(!file.exists(file.path(out, "Data", "RNA_09_CleanData", "RNA_CountData.RData"))){
  download.file("https://medical-epigenomics.org/papers/jakstat2024/Data.zip", destfile = file.path(out, "Data.zip"))
  unzip(file.path(out, "Data.zip"), exdir = out)
}

# load in vivo data at homeostasis
base::load(file.path(out, "Data", "RNA_09_CleanData", "RNA_CountData.RData"))
counts.in <- RNA.counts
gmap.in <- RNA.gmap
ann.in <- read_tsv(file.path(out, "Data", "RNA_09_CleanData", "Annotation.tsv"))

# load ex vivo data and subset to untreated
base::load(file.path(out, "Data", "TREAT_RNA_09_CleanData", "RNA_CountData.RData"))
counts.ex <- RNA.counts
gmap.ex <- RNA.gmap
ann.ex <- read_tsv(file.path(out, "Data", "TREAT_RNA_09_CleanData", "Annotation.tsv"))
ann.ex <- ann.ex %>%
  filter(treatment == "UT")
counts.ex <- counts.ex[,ann.ex$sample_name]


# combine data ------------------------------------------------------------

# combine matrices
stopifnot(all(row.names(counts.ex) == row.names(counts.in)))
gmap <- gmap.in
counts <- cbind(counts.ex, counts.in)

# combine annotations
colnames(ann.in)
col_names = c("sample_name", "experiment_id", "genotype", "cell_type", "treatment")
ann <- rbind(
  ann.in %>%select(any_of(col_names)),
  ann.ex %>% select(any_of(col_names))
) %>% 
  mutate(treatment = ifelse(treatment == "H", "in_vivo", "ex_vivo"))
ann <- data.table(ann)

write_rds(counts, file = file.path(out, "DEG_Counts.RDS"))
write_rds(gmap, file = file.path(out, "DEG_GMP.RDS"))
write_rds(ann, file = file.path(out, "DEG_Annotation.RDS"))



# DEG -----------------------------------------------------

res <- data.table()
(ctx <- ann$cell_type[2])
for(ctx in c("T8", "M")) {
  message("Analyzing ", ctx)
  annG <- ann[cell_type == ctx]
  
  # Count the number of replicates per genotype and treatment
  replicate_counts <- annG[, .N, by = c("genotype", "treatment")]
  
  # Identify genotypes with at least 3 replicates in both treatments (in_vivo and ex_vivo)
  valid_genotypes <- replicate_counts[replicate_counts$treatment == "in_vivo", 
                                      .(genotype, N_in_vivo = N)]
  valid_genotypes <- merge(valid_genotypes, 
                           replicate_counts[replicate_counts$treatment == "ex_vivo", 
                                            .(genotype, N_ex_vivo = N)], by = "genotype")
  
  # Filter genotypes with at least 3 replicates in both treatments
  valid_genotypes <- valid_genotypes[N_in_vivo >= 2 & N_ex_vivo >= 2]
  
  # If there are no valid genotypes, skip to the next cell type
  if(nrow(valid_genotypes) == 0) {
    next
  }
  
  # Keep only the valid genotypes
  annG <- annG[genotype %in% valid_genotypes$genotype]
  
  # Prepare annotation data
  annData <- data.frame(annG[, c("genotype", "experiment_id", "treatment"), with = F], row.names = annG$sample_name)
  annData$genotype <- factor(annData$genotype, levels = c("WT", unique(setdiff(annData$genotype, "WT"))))
  annData$treatment <- factor(annData$treatment, levels = c("in_vivo", "ex_vivo"))
  
  # Preprocess counts data
  datG <- counts[, annG$sample_name]
  if (NORMALIZE) datG <- calcNormFactors(DGEList(datG))
  threshold <- 10
  drop <- which(apply(cpm(datG), 1, max) < threshold)
  datG <- datG[-drop, ]
  
  # Define the model formula (including interaction and main effects)
  dream.form <- ~ genotype * treatment + (1|experiment_id)
  desMat <- model.matrix(~ genotype * treatment + experiment_id, data = annData)
  
  # clean names
  colnames(desMat) <- make.names(colnames(desMat))
  # Check which coefficients are not estimable
  non_estimable_coefs <-  nonEstimable(desMat)
  
  # Remove non-estimable coefficients from the model matrix
  
  estimable_coefs <- setdiff(colnames(desMat), non_estimable_coefs)
  desMat <- desMat[, estimable_coefs, drop = FALSE]
  # Perform Voom transformation
  cleanDev(); pdf(file.path(out, paste0("Voom_", ctx, "_Before.pdf")), w = 6, h = 6)
  voomRes <- voom(datG, design = desMat, plot = TRUE)
  dev.off()
  
  relevant.genes <- names(which(apply(voomRes$E, 1, mean) > CUTOFFS$DE.RNA$AveExpr))
  
  cleanDev(); pdf(file.path(out, paste0("Voom_", ctx, "_After.pdf")), w = 6, h = 6)
  voomRes <- voom(datG[relevant.genes,], design = desMat, plot = TRUE)
  dev.off()
  
  dreamRes <- voomWithDreamWeights(
    counts = datG[relevant.genes,], 
    formula = dream.form, 
    data = annData,
    quiet = TRUE
  )
  
    
  dream.weights    <- if(LMM.VoomWeights) voomRes$weights else dreamRes$weights
  dream.expression <- if(LMM.VoomWeights) as.matrix(voomRes) else as.matrix(dreamRes)
  summary(as.vector(dream.weights))
  # Model fitting
  fit <- dream(
    quiet = TRUE, #BPPARAM = SerialParam(),
    exprObj = dream.expression, 
    formula = dream.form, 
    data = annData,
    REML = LMM.REML,
    useWeights = TRUE, 
    weightsMatrix = dream.weights)
  
  # Export design matrix
  write.csv(fit$design, file.path(out, paste0("Design_", ctx, ".csv")), quote = F, row.names = F)
  cleanDev(); pdf(file.path(out, paste0("Design_", ctx, ".pdf")), w = 10, h = 20)
  pheatmap(fit$design)
  dev.off()
  
  # Export results (interaction and genotype effects)
  coefs <- colnames(coef(fit))
  #coefs <- coefs[coefs != "(Intercept)"]
  
  # Extracting coefficients related to genotype effects and interaction
  # Look for genotype-related contrasts (KO-WT in different conditions)
  
  for (coefx in coefs) {
    res <- rbind(res, data.table(coef = gsub("", "", coefx), cell_type = ctx, 
                                 topTable(fit = fit, coef = coefx, number = nrow(counts)), keep.rownames = T))
  }
  
  unique(coefs)
  # Extract genotype-only coefficients (without interaction terms)
  genotype_coefs <- grep("^genotype", coefs, value = TRUE)
  
  # Extract interaction coefficients (containing ":treatmentex_vivo")
  interaction_coefs <- grep(":treatmentex_vivo$", coefs, value = TRUE)
  
  # Remove interaction terms from the genotype list
  genotype_coefs <- setdiff(genotype_coefs, interaction_coefs)
  ###################
  #
  #
  ct.mt <- data.table()
  col_names <- c()
  
  for (gen in genotype_coefs) {
    L1 <- getContrast(dreamRes, dream.form, annData, c(paste0(gen, ":treatmentex_vivo"), gen))
    L1[L1 == -1] <- 1
    
    col_names <- c(col_names, paste0("ex.vivo_", gen))
    ct.mt <- cbind(ct.mt, L1)  # Bind as columns
  }
  ct.mt <- as.matrix(ct.mt)
  colnames(ct.mt) <- col_names
  rownames(ct.mt) <- coefs
  plotContrasts(t(ct.mt)) 
  ggsave(file = file.path(out,"contrast.png"))
  #fitting contrasts
  fit.mt = dream(dreamRes, dream.form, annData, ct.mt
              )
  
  ## Initialize results table
  Res.contrast <- data.table()
  
  # Loop through each contrast
  for (contrast in colnames(ct.mt)) {
    
    # Extract results for the contrast
    top_res <- topTable(fit.mt, coef = contrast, number = Inf)  # Get all significant genes
    
    # Convert to data.table and add row names (gene names)
    top_res_dt <- data.table(rn = rownames(top_res), top_res)
    
    # Add contrast name as a column
    top_res_dt[, coef := contrast]
    top_res_dt[,cell_type := ctx,]
    # Append to results table
    Res.contrast <- rbind(Res.contrast, top_res_dt, fill = TRUE)
  }
  Res <- rbind(res,Res.contrast)
  write.tsv(Res, file.path(out, paste0("DEG_results", ctx,".tsv")))
}

write.tsv(Res, file.path(out, paste0("DEG_results", ".tsv")))
#########################################
############

# Initialize the result table
gsea.res <- data.table() 

run_gsea <- function(limmaRes, enr.terms, celltypes = NULL, coefs = NULL) {
  
  # Initialize the result table
  gsea_res <- data.table() 
  
  # Determine cell types to process
  if (is.null(celltypes)) {
    celltypes <- unique(limmaRes$cell_type)
  }
  
  # Determine coefficients to process
  if (is.null(coefs)) {
    coefs <- unique(limmaRes$coef)
  }
  
  # Loop through each cell type
  for (ct in celltypes) {
    
    # Loop through each coefficient
    for (de_grp in coefs) {
      
      # Loop through each database in the enrichment terms
      for (dbx in names(enr.terms)) {
        
        # Subset the limma results based on the current cell type and coefficient
        subset_limmaRes <- limmaRes[limmaRes$cell_type == ct & limmaRes$coef == de_grp, ]
        
        # Extract statistics (logFC) and assign gene names as names
        stats <- with(subset_limmaRes, setNames(logFC, nm = gene))
        
        # Skip this iteration if there are missing values in stats
        if (any(is.na(stats))) {
          next
        }
        
        # Perform fgsea analysis
        fgsea_output <- fgsea(
          pathways = enr.terms[[dbx]],
          stats = stats
          #minSize = 15,   # Example additional arguments, adjust as necessary
          #maxSize = 500,  # Example additional arguments, adjust as necessary
          #nperm = 1000    # Example additional arguments, adjust as necessary
        )
        
        # Check if fgsea output is not empty and append the results to gsea_res
        if (length(fgsea_output) > 0) {
          gsea_res <- rbind(gsea_res, data.table(fgsea_output,
                                                 coef = de_grp,
                                                 celltype = ct,
                                                 db = dbx))
        }
      }
    }
  }
  
  # Return the combined GSEA results
  return(gsea_res)
}

databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
enr.terms <- enrichrGetGenesets(databases)
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
gsea.res <- run_gsea(res, enr.terms, celltypes = unique(res$cell_type),
                     coefs =unique(res$coef))
gsea.res %>% write_rds(out("fgsea_hom_vs_ex.vivo_per_CT.rds"))