require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)

# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/"
dir.create(out)


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
load(file.path(out, "Data", "RNA_09_CleanData", "RNA_CountData.RData"))
counts.in <- RNA.counts
gmap.in <- RNA.gmap
ann.in <- read_tsv(file.path(out, "Data", "RNA_09_CleanData", "Annotation.tsv"))

# load ex vivo data and subset to untreated
load(file.path(out, "Data", "TREAT_RNA_09_CleanData", "RNA_CountData.RData"))
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


# DEG ---------------------------------------------------------------------

res <- data.table()
(ctx <- ann$cell_type[1])
for(ctx in unique(ann$cell_type)){
  message("Analyzing ", ctx)
  annG <- ann[cell_type == ctx]
  # filter out genotypes, where there is not enough untreated / wild-type data
  annG.gts <- with(annG[genotype != "WT"], split(experiment_id, genotype))
  # gtx <- "STAT1KO"
  (genotypes.use <- names(which(sapply(names(annG.gts), function(gtx){
    x <- annG.gts[[gtx]]
    nrow(annG[(genotype == "WT" | genotype == gtx) & experiment_id %in% x][, .N, by=c("genotype", "treatment")][N >= 1])
  }) == 4)))
  annG <- annG[genotype %in% c("WT", genotypes.use)]
  
  annData <- data.frame(annG[,c("genotype", "experiment_id", "treatment"),with=F], row.names=annG$sample_name)
  annData$genotype <- factor(annData$genotype, levels=c("WT", unique(setdiff(annData$genotype, "WT"))))
  annData$treatment <- factor(annData$treatment, levels=c("in_vivo", "ex_vivo"))
  
  ggplot(data.table(annData)[,.N,by=c("genotype", "experiment_id", "treatment")], 
         aes(x=experiment_id, y=genotype, fill=N)) +
    facet_wrap(~treatment) + 
    geom_tile() + xRot()
  ggsave(file.path(out, paste0("CrossTable_", ctx, ".pdf")), w=12,h=6)
  
  datG <- counts[,annG$sample_name]
  if(NORMALIZE) datG <- calcNormFactors(DGEList(datG))
  dream.form <- ~ genotype * treatment + (1|experiment_id)
  desMat <- model.matrix(~ genotype*treatment + experiment_id, data=annData)
  
  cleanDev(); pdf(file.path(out, paste0("Voom_", ctx, "_Before.pdf")), w=6,h=6)
  voomRes <- voom(datG, design = desMat,plot=TRUE)
  dev.off()
  
  relevant.genes <- names(which(apply(voomRes$E, 1, mean) > CUTOFFS$DE.RNA$AveExpr))
  
  cleanDev(); pdf(file.path(out, paste0("Voom_", ctx, "_After.pdf")), w=6,h=6)
  voomRes <- voom(datG[relevant.genes,], design = desMat,plot=TRUE)
  dev.off()
  
  dreamRes <- voomWithDreamWeights(
    counts = datG[relevant.genes,], 
    formula = dream.form, 
    data = annData,
    quiet=TRUE, 
    BPPARAM=BPPARAM)
  
  dream.weights    <- if(LMM.VoomWeights) voomRes$weights else dreamRes$weights
  dream.expression <- if(LMM.VoomWeights) as.matrix(voomRes) else as.matrix(dreamRes)
  
  # model fit
  fit <- dream(
    quiet=TRUE, BPPARAM=BPPARAM,
    exprObj=dream.expression, 
    formula=dream.form, 
    data=annData,
    REML=LMM.REML,
    useWeights=TRUE, 
    weightsMatrix=dream.weights)
  
  # plot and write design
  write.csv(fit$design, file.path(out, paste0("Design_", ctx, ".csv")), quote=F, row.names=F)
  cleanDev(); pdf(file.path(out, paste0("Design_", ctx, ".pdf")), w=10,h=20)
  pheatmap(fit$design)
  dev.off()
  
  # Export results
  coefs <- colnames(coef(fit))
  coefs <- coefs[coefs != "(Intercept)"]
  coefx <- coefs[1]
  for(coefx in grep("", coefs, value=T)){
    res <- rbind(res, data.table(genotype=gsub("", "", coefx), cell_type = ctx, topTable(fit=fit,coef=coefx, number=nrow(counts)), keep.rownames=T))
  }
}
res[, B := NA]
res$z.std <- NULL
write.tsv(res, file.path(out, paste0("DEG", ".tsv")))


# Visualize results ---------------------------------------------------------
ggplot(res, aes(x=P.Value)) + geom_histogram(bins=20) + 
  facet_wrap(genotype ~ cell_type, scale="free_y")
ggsave(file.path(out, "DEG_Pvalue_distribution.pdf"), w=25,h=25)
  
ggplot(res, aes(x=logFC, y=-log10(P.Value))) + stat_binhex(aes(fill=log10(..count..))) + 
  facet_wrap(genotype ~ cell_type, scale="free_y")
ggsave(file.path(out, "DEG_Vulcanos.pdf"), w=10,h=10) 

with(res[adj.P.Val < 0.05], table(genotype, cell_type))

ggplot(res[adj.P.Val < 0.05][genotype != "treatmentex_vivo"], aes(x=genotype)) + 
  geom_bar() + xRot() +
  facet_grid(cell_type ~ .) +
  scale_y_log10()
ggsave(file.path(out, "DEG_Numbers.pdf"), w=10,h=10)