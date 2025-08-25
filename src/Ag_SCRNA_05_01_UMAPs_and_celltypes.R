source("src/00_init.R")
out <- dirout("Ag_SCRNA_05_01_UMAPs_and_celltypes")
require(pdist)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")
registerDoMC(cores=10)
InDir <- dirout("Ag_SCRNA_02_01_Integration")
InDir1 <- dirout("Ag_SCRNA_04_01_proj_ex.vivo/")
# Original UMAPs ----------------------------------------------------------
mobjs <- list()
for(tissuex in tissue <-c("ex.vivo","in.vivo")){
  (base::load(InDir(paste0("in.vivo","/soupx/MonocleObject.RData"))))
  mobjs[[tissuex]] <- monocle.obj
}


# Duplet scors ------------------------------------------------------------


ff <- list.files("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/Ag_SCRNA_03_01_Duplets", 
                             pattern="Duplet_Scores_.*.tsv", full.names = TRUE)
ds <- lapply(ff, fread)
names(ds) <- basename(ff)
ds <- rbindlist(ds, idcol = "sample")
ds[, sample := gsub("Duplet_Scores_(.+).tsv", "\\1", sample)]
ds[, rn := paste0(rn, "_", sample)]
ds <- ds[,c("rn", "dupletScore")]


# Projection from Monocle (original) --------------------------------------
res <- data.table()
tx <- names(mobjs)[1]
for(tx in names(mobjs)){
  monocle.obj <- mobjs[[tx]]
  dDT.umap <- cbind(
    data.table(reducedDims(monocle.obj)$UMAP, keep.rownames = TRUE),
    data.table(
      sample=monocle.obj$sample,
      tissue=monocle.obj$tissue
    ))
  res <- rbind(res, dDT.umap)
}
colnames(res)[colnames(res) %in% c("V1", "V2")] <- c("UMAP_1", "UMAP_2")
stopifnot(!any(is.na(res$UMAP_1)))

# Duplets
res <- merge(res, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(res[is.na(dupletScore)]) == 0)
table(res$dupletScore > 0.9)
ggplot(res[!is.na(tissue)], aes(x=UMAP_1, y=UMAP_2)) + 
  scale_fill_hexbin() +
  facet_wrap(~tissue, ncol=3) +
  theme_bw(12) +
  stat_summary_hex(aes(z=dupletScore),fun=mean, bins=100)
ggsave(out("ProjMonocle_Duplets.pdf"), w=12,h=4)
table(res$tissue)
#res <- res[dupletScore < 0.9 | tissue != "in.vivo"]
res <- res[dupletScore < 0.9]

# Save
saveRDS(res, out("ProjMonocle.RDS"))

# Clusters
mnam <- names(mobjs)[1]
dDT.ct <- list()
for(mnam in names(mobjs)){
  monocle.obj <- mobjs[[mnam]]
  dDT.ct[[mnam]] <- data.table(data.frame(colData(monocle.obj)@listData), keep.rownames = TRUE)[,c("rn", "sample")][match(colnames(monocle.obj), rn)]
  dDT.ct[[mnam]]$functional.cluster <- getCL(monocle.obj)
  dDT.ct[[mnam]]$functional.cluster.conf <- 1
}
dDT.ct <- rbindlist(dDT.ct)
dDT.ct <- dDT.ct[match(res$rn, rn)]
stopifnot(!any(is.na(dDT.ct$sample)))
saveRDS(dDT.ct, out("ProjMonocle_Clusters.RDS"))




# celltypes from singleR after further curation -------------------------------------------------------


# Projection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("Ag_SCRNA_proj_ex.vivo_inc_Mye/")(""), pattern="Output_ex", full.names = TRUE)
dL <- lapply(ff, fread)

dDT <- rbindlist(dL)
head(dDT)
invivo <- fread(InDir1("Output_in.vivo.tsv"))
invivo$rn <- invivo$cell_id
dDT <- rbind(dDT, invivo, fill = T)

dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
#
dDT.umap <- merge(dDT.umap, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(dDT.umap[is.na(dupletScore)]) == 0)
table(is.na(dDT.umap$rn))

#dDT.umap <- dDT.umap[dupletScore < 0.9 | tissue != "in.vivo"]
#dDT.umap <- dDT.umap[dupletScore < 0.9]
table(dDT.umap$tissue)
saveRDS(dDT.umap, out("ProjVivo.RDS"))

dDT.ct <- dDT[,c("rn", "sample", "functional.cluster", "functional.cluster.conf"), with=F]
dDT.ct <- dDT.ct[match(dDT.umap$rn, rn)]
stopifnot(!any(is.na(dDT.ct$sample)))
saveRDS(dDT.ct, out("ProjVivo_celltypes.RDS"))

# Crossprojection_Invivo -------------------------------------------------------
ff <- list.files(dirout_load("Ag_SCRNA_proj_ex.vivo_inc_Mye/")(""), pattern="OutputCrossprojection_", full.names = TRUE)
dL <- lapply(ff, fread)
dDT <- rbindlist(dL, fill=TRUE)
dDT.umap <- dDT[,c("rn", "sample", "UMAP_1", "UMAP_2", "tissue"), with=F]
dDT.umap.invivo <- readRDS(out("ProjMonocle.RDS"))[,-"dupletScore", with=F]
dDT.umap <- rbind(dDT.umap, dDT.umap.invivo[tissue == "in.vivo"])
dDT.umap <- merge(dDT.umap, ds, by="rn", all=TRUE)[!is.na(UMAP_1)]
stopifnot(nrow(dDT.umap[is.na(dupletScore)]) == 0)
table(dDT.umap$tissue)
#dDT.umap <- dDT.umap[dupletScore < 0.9 | tissue != "in.vivo"]
dDT.umap <- dDT.umap[dupletScore < 0.9]
table(dDT.umap$tissue)
saveRDS(dDT.umap, out("ProjVivoX.RDS"))

