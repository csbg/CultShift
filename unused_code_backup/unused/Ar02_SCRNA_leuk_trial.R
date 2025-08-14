source("src/00_init.R")
require(ProjecTILs)
require(umap)
require(ggrepel)


out1 <- dirout("SCRNA_50_01_Trajectories/")
source("src/FUNC_ProjecTILs_PLUS_Ar1.R")
base.dir <- "SCRNA_50_02_ProjectionTrajectories/"
out <- dirout(base.dir)
SANN <- fread(PATHS$SCRNA$ANN)
ann1 <- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle.singleR")("Annotation.tsv"))
ann <- fread(dirout_load("SCRNA_20_Summary/in.vivo_monocle.singleR")("Annotation.tsv"))

basedir_corr <- "SCRNA_33_DE_Nebula_testClustering/Output/exvivo_9dinvivo14d/"
out_corr <- dirout(basedir_corr)
correlation_logfc<-fread(out_corr("correlation_logfc_invivo_14dvs_leukvivo9d.tsv"))
colnames(correlation_logfc)<-c("gene","correlation")

# test statistcs ---------------------------------------------------------------
ff <- list.files(out1(""), pattern="Pseudotime.*.tsv")
names(ff) <- gsub("Pseudotime_(.+).tsv", "\\1", ff)

pDT <- rbindlist(lapply(ff, function(fx) fread(out1(fx))), idcol = "celltype")
pDT[, traj := pseudotime.cds.]
pDT$pseudotime.cds. <- NULL
write.tsv(pDT, out1("Values_in.vivo.tsv"))
pDT <- merge(pDT, ann[Clusters != "HSC"], by="rn")[!is.na(mixscape_class.global)]
pDT[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT[, traj.scale := scale(traj), by="celltype"]
pDT <- pDT[timepoint != "28d"]

# . test_in.vivo ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res <- data.table()
for(typex in unique(pDT$celltype)){
  pDT1 <- pDT[celltype == typex]
  for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT1[gene == gx]$traj
    x2 <- pDT1[gene == "NTC"]$traj
    if(length(x1) > 10 & length(x2) > 10){
      res <- rbind(res, data.table(
        p.wx=wilcox.test(x1, x2)$p.value,
        p.ks=ks.test(x1, x2)$p.value,
        d=median(x1) - median(x2),
        type=typex,
        gene=gx
      ))
    }
  }
}
res[, padj.wx := p.adjust(p.wx, method="BH")]
res[, padj.ks := p.adjust(p.ks, method="BH")]
write.tsv(res, out1("Statistics_in.vivo.tsv"))


# . load ------------------------------------------------------------------
res <- fread(out1("Statistics_in.vivo.tsv"))


# . plot stats_in.vivo ------------------------------------------------------------
ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out1("Statistics_Comparison_in.vivo.pdf"),w=20,h=6)

# xDT <- melt(res, id.vars = c("type", "gene"))
# xDT[, measurement := gsub("\\..+$", "", variable)]
# xDT[, type := gsub("^.+?\\.", "", variable)]
ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out1("Statistics_in.vivo.pdf"), w=4,h=10)



# . UMAP ------------------------------------------------------------------
pDT.UMAP <- pDT[abs(traj.scale) < 3]
ggplot(pDT.UMAP, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out1("UMAP_in.vivo.pdf"), w=5,h=5)


# . plot distributions ----------------------------------------------------
pDT.distr <- copy(pDT)

pDT.distr<-pDT.distr[gene %in% c(unique(res$gene),"NTC")]
ggplot(pDT.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
pDT.distr <- pDT.distr[abs(traj.scale) < 3]
pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("tissue","gene", "celltype")]
pDT.sum[, traj_diff := traj]

cellt<-"Ery"
gen="j"
#Subtract NTC from all gene and store it i traj_diff column
for (cellt in unique(pDT.sum$celltype)){
  for (gen in unique(pDT.sum$gene)){
    pDT.sum[celltype==cellt & gene == gen]$traj_diff = pDT.sum[celltype==cellt & gene == gen]$traj_diff-pDT.sum[celltype == cellt & gene == "NTC"]$traj}
}


pDT.stats <- copy(res)
pDT.stats[, celltype := type]
pDT.stats[, type := "not.sig"]
pDT.stats[padj.ks <= 0.1, type := "sig.low"]
pDT.stats[padj.ks <= 0.01, type := "sig.high"]
pDT.distr <- merge(pDT.distr, pDT.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)

pDT.distr[gene=="NTC", type := "NTC"]
pDT.distr<-pDT.distr[gene %in% c(pDT.stats$gene,"NTC")]


ggplot(pDT.distr, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum, color="black") + 
  geom_errorbarh(data=pDT.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(out1("Distribution_in.vivo.pdf"), w=20,h=10)

ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggsave(out1("Distribution_Brd9.pdf"), w=5,h=4)


# . scatterplot Ery_vs_Mye in.vivo
-----------------------------------------------------------
pDT2 <- pDT[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
pDT2[, V1 := scale(V1), by="celltype"]
pDT2 <- dcast.data.table(pDT2, gene ~ celltype, value.var = "V1")
ggplot(pDT2, aes(x=Ery, y=Mye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out1("Scatter_EryVsMye_in.vivo.pdf"),w=8,h=8)




# Annotation --------------------------------------------------------------



ff3 <- list.files(out(""), pattern="Output.*.tsv")
ff3 <- ff3[!grepl("in.vivo", ff3)]
ff3 <- ff3[!grepl("ex.vivo", ff3)]
names(ff3) <- gsub("Output_(.+).tsv", "\\1", ff3)
ff3
pDT_leuk <- rbindlist(lapply(ff3, function(fx) fread(out(fx))), idcol = "sample")
pDT_leuk[, traj := pseudotime]
pDT_leuk$pseudotime <- NULL

#write.tsv(pDT_leuk, out("Values_leuk.vivo.tsv"))
pDT_leuk <- merge(pDT_leuk, ann1, by="rn")[!is.na(mixscape_class.global)]
pDT_leuk[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT_leuk[, celltype := ct]
pDT_leuk[, traj.scale := scale(traj), by="celltype"]
pDT_leuk[, id := paste0(tissue, "_", timepoint)]
unique(pDT_leuk$id)
(idx <- unique(pDT_leuk$id)[1])
for(idx in unique(pDT_leuk$id)){
  outS <- dirout(paste0(base.dir, idx))
  pDT_leuk <- pDT_leuk[id == idx]
  
  # . test ------------------------------------------------------------------
  typex <- "Ery"
  gx <- "Rcor1"
  res <- data.table()
  for(typex in unique(pDT_leuk$celltype)){
    pDT_leuk1 <- pDT_leuk[celltype == typex]
    for(gx in unique(pDT_leuk1[mixscape_class.global != "NTC"]$gene)){
      x1 <- pDT_leuk1[gene == gx]$traj
      x2 <- pDT_leuk1[gene == "NTC"]$traj
      if(length(x1) > 10 & length(x2) > 10){
        res <- rbind(res, data.table(
          p.wx=wilcox.test(x1, x2)$p.value,
          p.ks=ks.test(x1, x2)$p.value,
          d=median(x1) - median(x2),
          type=typex,
          gene=gx
        ))
      }
    }
  }
  res[, padj.wx := p.adjust(p.wx, method="BH")]
  res[, padj.ks := p.adjust(p.ks, method="BH")]
  write.tsv(res, outS("Statistics.tsv"))
  
  
  # . load ------------------------------------------------------------------
  #res <- fread(outS("Statistics.tsv"))
  
  # . plot stats ------------------------------------------------------------
  ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
    geom_point() +
    geom_point(shape=1, color="black") +
    theme_bw(12)+ 
    geom_text_repel(aes(label=paste(gene)), color="black")+
    scale_color_gradient2(low="blue", high="red") +
    facet_grid(. ~ type)
  ggsave(outS("Statistics_Comparison.pdf"),w=20,h=6)
  
  # xDT <- melt(res, id.vars = c("type", "gene"))
  # xDT[, measurement := gsub("\\..+$", "", variable)]
  # xDT[, type := gsub("^.+?\\.", "", variable)]
  ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
    geom_point() +
    theme_bw(12)+ 
    scale_color_gradient2(low="blue", high="red") +
    xRot()
  ggsave(outS("Statistics.pdf"), w=4,h=10)
  
  # . UMAP ------------------------------------------------------------------
  pDT_leuk.UMAP <- pDT_leuk[abs(traj.scale) < 3]
  ggplot(pDT_leuk.UMAP, aes(x=UMAP1, y=UMAP2)) + 
    stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
    theme_bw(12) +
    scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
  ggsave(outS("UMAP.pdf"), w=5,h=5)
  
  
  # . celltypes -------------------------------------------------------------
  ggplot(pDT_leuk, aes(x=UMAP1, y=UMAP2, color=celltype)) + 
    geom_point() +
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP.jpg"), w=5,h=5)
  
  ggplot(pDT_leuk, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex() +
    facet_grid(. ~ celltype) + 
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP_hex.pdf"), w=20,h=5)
  
  ggplot(pDT_leuk, aes(x=celltype)) + 
    geom_bar() +
    scale_y_log10() +
    theme_bw(12)
  ggsave(outS("Celltypes_Nubmers.pdf"), w=5,h=5)
  
  # . plot distributions ----------------------------------------------------
  pDT_leuk.distr <- copy(pDT_leuk)
  ggplot(pDT_leuk.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
  pDT_leuk.distr <- pDT_leuk.distr[abs(traj.scale) < 3]
  pDT_leuk.sum <- pDT_leuk.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene", "celltype")]
  pDT_leuk.stats <- copy(res)
  pDT_leuk.stats[, celltype := type]
  pDT_leuk.stats[, type := "not.sig"]
  pDT_leuk.stats[p.ks < 0.1, type := "sig.low"]
  pDT_leuk.stats[p.ks < 0.01, type := "sig.high"]
  pDT_leuk.distr <- merge(pDT_leuk.distr, pDT_leuk.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)
  pDT_leuk.distr[gene == "NTC", type := "NTC"]
  ggplot(pDT_leuk.distr, aes(y=gene, x=traj)) + 
    geom_violin(color=NA, aes(fill=type), scale="width") + 
    scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
    geom_point(data=pDT_leuk.sum, color="black") + 
    geom_errorbarh(data=pDT_leuk.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
    theme_bw(12) + 
    facet_grid(. ~ celltype, scale="free_x") +
    xRot()
  ggsave(outS("Distribution_leuk.vivo.pdf"), w=20,h=10)
  
  ggplot(pDT_leuk.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
         aes(x=traj, color=gene)) + 
    theme_bw(12) +
    geom_density()
  ggsave(outS("Distribution_Brd9_leuk.vivo.pdf"), w=5,h=4)
  
  
  # . scatterplot -----------------------------------------------------------
  pDT_leuk2 <- pDT_leuk[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
  pDT_leuk2[, V1 := scale(V1), by="celltype"]
  pDT_leuk2 <- dcast.data.table(pDT_leuk2, gene ~ celltype, value.var = "V1")
  ggplot(pDT_leuk2, aes(x=Ery, y=Mye)) + 
    geom_point() + 
    theme_bw(12) +
    geom_abline() + 
    geom_text_repel(aes(label=gene))
  ggsave(outS("Scatter_EryVsMye_leuk.vivo.pdf"),w=8,h=8)
}
#test leukemia-------------------------------------------
# . test ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res_leuk <- data.table()
for(typex in unique(pDT_leuk$celltype)){
  pDT1 <- pDT_leuk[celltype == typex]
  for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT1[gene == gx]$traj
    x2 <- pDT1[gene == "NTC"]$traj
    if(length(x1) > 10 & length(x2) > 10){
      res_leuk <- rbind(res_leuk, data.table(
        p.wx=wilcox.test(x1, x2)$p.value,
        p.ks=ks.test(x1, x2)$p.value,
        d=median(x1) - median(x2),
        type=typex,
        gene=gx
      ))
    }
  }
}
res_leuk[, padj.wx := p.adjust(p.wx, method="BH")]
res_leuk[, padj.ks := p.adjust(p.ks, method="BH")]
write.tsv(res_leuk, outSL("Statistics_leuk.vivo.tsv"))
res_leuk <- fread(outSL("Statistics_leuk.vivo.tsv"))
ggplot(res_leuk, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison_leuk.vivo.pdf"),w=20,h=6)

#
ggplot(res_leuk, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_leuk.vivo.pdf"), w=4,h=10)

# . UMAP ------------------------------------------------------------------
pDT.UMAP_leuk <- pDT_leuk[abs(traj.scale) < 3]
ggplot(pDT.UMAP_leuk, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_leuk.vivo.pdf"), w=5,h=5)
#

# . plot distributions leukemia----------------------------------------------------
pDT.distr_leuk <- copy(pDT_leuk)
ggplot(pDT.distr_leuk, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
pDT.distr_leuk <- pDT.distr_leuk[abs(traj.scale) < 3]
pDT.sum_leuk <- pDT.distr_leuk[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("tissue","gene", "celltype")]

pDT.sum_leuk[, traj_diff := traj]

cellt<-"Ery"
gen="j"
for (cellt in unique(pDT.sum_leuk$celltype)){
  for (gen in unique(pDT.sum_leuk$gene)){
    pDT.sum_leuk[celltype==cellt & gene == gen]$traj_diff = pDT.sum_leuk[celltype==cellt & gene == gen]$traj_diff-pDT.sum_leuk[celltype == cellt & gene == "NTC"]$traj}
}

pDT.stats_leuk <- copy(res_leuk)
pDT.stats_leuk[, celltype := type]
pDT.stats_leuk[, type := "not.sig"]
pDT.stats_leuk[padj.ks <= 0.1, type := "sig.low"]
pDT.stats_leuk[padj.ks <= 0.01, type := "sig.high"]
pDT.distr_leuk <- merge(pDT.distr_leuk, pDT.stats_leuk[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)

pDT.distr_leuk[gene=="NTC", type := "NTC"]

ggplot(pDT.distr_leuk[gene=="NTC"&celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()

ggplot(pDT.distr_leuk, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_leuk, color="black") + 
  geom_errorbarh(data=pDT.sum_leuk, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(outS("Distribution_leuk.vivo.pdf"), w=20,h=10)
########################################################
ggplot(pDT_ex.distr, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum[gene %in% pDT.distr_leuk$gene & celltype =="Mye"], color="black") + 
  geom_errorbarh(data=pDT.sum[gene %in% pDT.distr_leuk$gene & celltype =="Mye"], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + xlim(0,40)+
  #facet_grid(. ~ celltype) +
  xRot()
ggsave(outS("Distribution_invivo_Mye.pdf"), w=20,h=10)


ggplot(pDT.distr_leuk[celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_leuk[gene %in% pDT.distr_leuk$gene & celltype =="Mye"], color="black") + 
  geom_errorbarh(data=pDT.sum_leuk[gene %in% pDT.distr_leuk$gene & celltype =="Mye"], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) +  xlim(0,40)+
  facet_grid(. ~ celltype) +
  xRot()
ggsave(outS("Distribution_leukvivo_Mye.pdf"), w=20,h=10)


############################################
#combine leukemia and in.vivo-------------------------

pDT.sum_leuk[gene=="NTC"]
pDT.sum_leuk_combined<-rbind(pDT.sum_ex[gene %in% pDT.distr_leuk$gene],pDT.sum_leuk)
pDT.sum_leuk_combined_Mye<-pDT.sum_leuk_combined[celltype=="Mye"]

Mye_leuk<-pDT.distr_leuk[celltype=="Mye"][,c("gene","rn","UMAP1","UMAP2","traj","sample.x","celltype","tissue","timepoint","type","traj.scale")]
colnames(Mye_leuk)<-gsub("sample.x","sample",colnames(Mye_leuk))

pDT_ex.distr<-pDT_ex.distr[,c("gene","rn","UMAP1","UMAP2","traj","sample","celltype","tissue","timepoint","type","traj.scale")]
pDT_ex.distr<-pDT_ex.distr[gene %in% pDT.distr_leuk$gene & celltype =="Mye"]
Mye_leuk<-rbind(Mye_leuk,pDT_ex.distr)


#combined_plts-------------------------------




ggplot(Mye_leuk[celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_leuk_combined_Mye[gene %in% pDT.distr_leuk$gene], color="black") + 
  geom_errorbarh(data=pDT.sum_leuk_combined_Mye[gene %in% pDT.distr_leuk$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_leuk.ex.pdf"), w=20,h=10)


ggplot(data=pDT.sum_leuk_combined_Mye[gene %in% pDT.distr_leuk$gene], aes(y=gene, x=traj_diff)) + 
  geom_point(data=pDT.sum_leuk_combined_Mye[gene %in% pDT.distr_leuk$gene], color="black") + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_leuk.ex_traj-ntc.pdf"), w=20,h=10)

Mye_leuk<-Mye_leuk[celltype=="Mye"]

#scatter not scaled

pDT2_com <- Mye_leuk[tissue %in% c("leukemia", "ex.vivo")][, median(traj), by=c("tissue", "gene")]

#pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")

ggplot(pDT2_com, aes(x=leukemia, y=ex.vivo)) + 
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(outS("Scatter_leuk_ex.vivo.pdf"),w=8,h=8)

#scatter scaled
pDT2_com <- Mye[tissue %in% c("leukemia", "ex.vivo")][, median(traj), by=c("tissue", "gene")]
pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")
head(pDT2_com)
ggplot(pDT2_com, aes(x=leukemia, y=ex.vivo)) + 
  geom_point() + 
  geom_abline()+
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(outS("Scatter_leuk_in.vivo_scaled.png"))



#traj_median(gene)-ntc

pDT_summary_leuk_Mye<-pDT.sum_leuk_combined_Mye[,c("tissue","gene","traj_diff")]

pDT_summary_leuk_Mye<- dcast.data.table(pDT_summary_leuk_Mye, gene ~ tissue, value.var = "traj_diff")
pDT_summary_leuk_Mye[,diff:=leukemia-in.vivo]
head(pDT_summary_leuk_Mye)
ggplot(pDT_summary_leuk_Mye, aes(x=leukemia, y=in.vivo)) + 
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_leuk_in.vivo_traj-ntc.pdf"),w=8,h=8)






#scale to sd(NTC)

Mye[, traj_ntc_norm := traj]
tis="leukemia"
gen="j"
for (tis in unique(Mye$tissue)){
  Mye[tissue==tis]$traj_ntc_norm = (Mye[tissue==tis]$traj-(mean(Mye[tissue==tis & gene == "NTC"]$traj)))/sd(Mye[tissue==tis & gene == "NTC"]$traj)
}



###########################
#summary-median of (gene_traj-mean_ntc)/sd(ntc)
pDT.ntc_leuk_sum<-  Mye[, .(median_ntc_norm_traj = median(traj_ntc_norm), q1 = quantile(traj_ntc_norm, 0.25), q2 = quantile(traj_ntc_norm, 0.75)),by=c("tissue","gene")]
pDT.ntc_leuk_sum[, traj_ntc_norm:= median_ntc_norm_traj]


pDT.ntc_leuk_sum_mean<-  Mye[, .(mean_ntc_norm_traj = mean(traj_ntc_norm), q1 = mean(traj_ntc_norm)-sd(traj_ntc_norm), q2 = mean(traj_ntc_norm)+sd(traj_ntc_norm)),by=c("tissue","gene")]
pDT.ntc_leuk_sum_mean[, traj_ntc_norm:= mean_ntc_norm_traj]




#median
ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_leuk_sum[gene %in% pDT.distr_leuk$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_leuk_sum[gene %in% pDT.distr_leuk$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(out("Distribution_Median_Mye_leuk_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)


ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_leuk_sum[gene %in% pDT.distr_leuk$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_leuk_sum[gene %in% pDT.distr_leuk$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()

head(pDT.ntc_leuk_sum)
pDT.ntc_leuk_sum_med<-pDT.ntc_leuk_sum[,c("tissue","gene","median_ntc_norm_traj")]
colnames(pDT.ntc_leuk_sum_med)<-c("tissue","gene","traj_ntc_norm")
head(pDT.ntc_leuk_sum_med)

pDT.ntc_leuk_sum_med_diff<- dcast.data.table(pDT.ntc_leuk_sum_med, gene ~ tissue, value.var = "traj_ntc_norm")
pDT.ntc_leuk_sum_med_diff<-pDT.ntc_leuk_sum_med_diff[gene!="Chaf1a"]
pDT.ntc_leuk_sum_med_diff[, diff:= in.vivo-leukemia]
head(pDT.ntc_leuk_sum_med_diff)

ggplot() + 
  #geom_segment(data=pDT.ntc_leuk_sum_med_diff, mapping=aes(x=leukemia, y=in.vivo, xend=leukemia+diff, yend=in.vivo), arrow = arrow(length = unit(0.05, "cm")), size=0.05, color="blue") + 
  geom_point(data=pDT.ntc_leuk_sum_med_diff, mapping=aes(x=leukemia, y=in.vivo), size=1, shape=21, fill="white") +
  theme_bw(12) +geom_text_repel(aes(label=gene))
ggsave(out("segment.pdf"),w=8,h=8)
head(pDT.ntc_leuk_sum_med_diff)

ggplot(pDT.ntc_leuk_sum_med_diff, aes(x=leukemia, y=gene)) + 
  geom_point(size=1, shape=21, fill="white")+
  geom_segment(data=pDT.ntc_leuk_sum_med_diff, mapping=aes(x=leukemia, y=gene, 
                                                      xend=leukemia+diff, yend=gene),arrow = arrow(length = unit(0.09, "cm")), 
               size=0.09, color="black")+
  theme_bw(12)
#geom_text_repel(aes(label=gene,max.overlaps =40))

ggsave(out("thenga.png"))
#mean
ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_leuk_sum_mean[gene %in% pDT.distr_leuk$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_leuk_sum_mean[gene %in% pDT.distr_leuk$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(out("Distribution_Mean_Mye_leuk_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)

##
pDT.ntc_leuk_sum_mean<-pDT.ntc_leuk_sum_mean[,c("tissue","gene","traj_ntc_norm")]
head(pDT.ntc_leuk_sum_mean)

pDT.ntc_leuk_sum_mean_diff<- dcast.data.table(pDT.ntc_leuk_sum_mean, gene ~ tissue, value.var = "traj_ntc_norm")
pDT.ntc_leuk_sum_mean_diff<-pDT.ntc_leuk_sum_mean_diff[gene!="Chaf1a"]
pDT.ntc_leuk_sum_mean_diff[, diff:= in.vivo-leukemia]
head(pDT.ntc_leuk_sum_mean_diff)

ggplot(pDT.ntc_leuk_sum_med_diff, aes(x=leukemia, y=in.vivo)) + 
  geom_point(size=2, shape=21, fill="white") +
  geom_segment(data=pDT.ntc_leuk_sum_mean_diff, mapping=aes(x=leukemia, 
                                                       y=in.vivo, xend=leukemia+diff, yend=in.vivo), 
               arrow = arrow(length = unit(0.2, "cm")), size=0.2, color="black") + 
  theme_bw(12)+ggrepel
ggsave(out("segment.pdf"),w=8,h=8)
########################################
#
pDT.ntc_leuk_sum<-pDT.ntc_leuk_sum[,c("gene","tissue","median_ntc_norm_traj")]

pDT.ntc_leuk_sum_plot<-dcast.data.table(pDT.ntc_leuk_sum, gene ~ tissue, value.var = "median_ntc_norm_traj")

pDT.ntc_leuk_sum_plot_color<-merge(pDT.ntc_leuk_sum_plot,correlation_logfc,by="gene")

NTC<-pDT.ntc_leuk_sum_plot[gene=="NTC"]
correlation<-0.00
NTC<-cbind(NTC,correlation)

pDT.ntc_leuk_sum_plot_color<-rbind(pDT.ntc_leuk_sum_plot_color,NTC)
tail(pDT.ntc_leuk_sum_plot_color)
#scatter_plot
pDT.ntc_leuk_sum_plot[gene=="NTC"]
pDTA<-pDT.ntc_leuk_sum_plot_color


#normal-with NTC
ggplot(pDT.ntc_leuk_sum_plot, aes(x=leukemia, y=in.vivo)) + 
  geom_point()+
  theme_bw(12) +
  geom_text_repel(aes(label=gene))

ggsave(out("Scatter_Mye_leuk_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=7,h=8)

#color-No-NTC
ggplot(pDTA, aes(x=leukemia, y=in.vivo)) + 
  geom_point(aes(color = pDTA$correlation))+
  
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_leuk_in.vivo_traj_ntc_mean_sd_norm_color_corrlogfc.pdf"),w=8,h=8)

#
