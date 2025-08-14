source("src/00_init.R")
require(ProjecTILs)
require(umap)
require(ggrepel)
require(monocle3)

source("src/FUNC_ProjecTILs_PLUS_Ar1.R")

InDir <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide/")
basedir <- "Ag_ScRNA_01_Trajectories_comparison_annotation_projection"
out <- dirout(basedir)

ann <- fread(InDir("invivo_Annotations.tsv"))
ann1 <- fread(InDir("exvivo_Annotations.tsv"))
# Load annotation ---------------------------------------------------------
monocle_obj <- read_rds(InDir("invivo_monocle_proj.rds"))
str(monocle_obj)

celltypes <- fread('metadata/FIGS_celltypes.tsv')
celltypes <- celltypes[Name %in% ann$celltype_projection]


typex <- "Mye"
for(typex in unique(celltypes[Type != "HSC"]$Type)){
  print(typex)
  cells <- ann[celltype_projection %in% celltypes[Type %in% c(typex, "HSC")]$Name]$rn
  
  cds <- monocle_obj[,cells]
  
  # recluster (not used but necessary)
  cds <- monocle3::cluster_cells(cds)
  
  # learn graph
  cds <- learn_graph(
    cds,
    verbose = TRUE,
    use_partition = FALSE,
    close_loop = TRUE)
  
  # order cells
  cds <- order_cells(cds, root_cells=ann[rn %in% colnames(cds)][celltype_projection == "HSC"]$rn)
  
  # Make plot
  plot_cells(cds,color_cells_by = "pseudotime")
  ggsave(out("Pseudotime_",typex,".jpg"), w=10,h=10)
  
  # export table
  write.tsv(data.table(data.frame(pseudotime(cds)), keep.rownames = TRUE), out("Pseudotime_",typex,".tsv"))
}

###########
#checked_until_here--1/27/25

# test statistcs ---------------------------------------------------------------
ff <- list.files(out(""), pattern="Pseudotime.*.tsv")
names(ff) <- gsub("Pseudotime_(.+).tsv", "\\1", ff)

pDT <- rbindlist(lapply(ff, function(fx) fread(out(fx))), idcol = "celltype")
pDT[, traj := pseudotime.cds.]
pDT$pseudotime.cds. <- NULL
write.tsv(pDT, out("Values_in.vivo.tsv"))
pDT <- merge(pDT, ann[Clusters != "HSC"], by="rn")[!is.na(mixscape_class.global)]
pDT[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT[, traj.scale := scale(traj), by="celltype"]
pDT <- pDT[timepoint != "28d"]

# . test ------------------------------------------------------------------
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
write.tsv(res, out("Statistics_in.vivo.tsv"))


# . load ------------------------------------------------------------------
res <- fread(out("Statistics_in.vivo.tsv"))


# . plot stats in.vivo------------------------------------------------------------
ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison_in.vivo.pdf"),w=20,h=6)


ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_in.vivo.pdf"), w=4,h=10)



# . UMAP ------------------------------------------------------------------
pDT.UMAP <- pDT[abs(traj.scale) < 3]
ggplot(pDT.UMAP, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_in.vivo.pdf"), w=5,h=5)


# . plot distributions ----------------------------------------------------
pDT.distr <- copy(pDT)
pDT.distr<-pDT.distr[gene %in% c(unique(res$gene),"NTC")]
ggplot(pDT.distr, aes(x=traj, y=traj.scale)) + 
  geom_hex() + 
  facet_wrap(~celltype, scales = "free")
pDT.distr <- pDT.distr[abs(traj.scale) < 3]
unique(pDT.sum$celltype)
pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("tissue","gene", "celltype")]
pDT.sum[, traj_diff := traj]

cellt<-"Ery"
gen="j"
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
ggsave(out("Distribution_in.vivo.pdf"), w=20,h=10)

ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggsave(out("Distribution_Brd9.pdf"), w=5,h=4)


# . scatterplot -----------------------------------------------------------
pDT2 <- pDT[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
pDT2[, V1 := scale(V1), by="celltype"]
pDT2 <- dcast.data.table(pDT2, gene ~ celltype, value.var = "V1")
ggplot(pDT2, aes(x=Ery, y=Mye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_EryVsMye_in.vivo.pdf"),w=8,h=8)
#############################################################################
##########
#end------
##################


#exvivo_trajectories----------------


# ex.vivo --------------------------------------------------------------
ff1 <- list.files(inp(""), pattern="Output.*.tsv")
ff1 <- ff1[!grepl("in.vivo", ff1)]
ff1 <- ff1[!grepl("leukemia", ff1)]
names(ff1) <- gsub("Output_(.+).tsv", "\\1", ff1)

pDT_ex <- rbindlist(lapply(ff1, function(fx) fread(out(fx))), idcol = "sample")
pDT_ex[, traj := pseudotime]
pDT_ex$pseudotime <- NULL

write.tsv(pDT_ex, out("Values_ex.vivo.tsv"))
pDT_ex <- merge(pDT_ex, ann1, by="rn")[!is.na(mixscape_class.global)]
pDT_ex[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT_ex[, celltype := ct]
pDT_ex[, traj.scale := scale(traj), by="celltype"]
pDT_ex[, id := paste0(tissue, "_", timepoint)]

(idx <- unique(pDT_ex$id)[1])
for(idx in unique(pDT_ex$id)){
  outS <- dirout(paste0(base.dir, idx))
  pDT_ex <- pDT_ex[id == idx]
  
  # . test ------------------------------------------------------------------
  typex <- "Ery"
  gx <- "Rcor1"
  res <- data.table()
  for(typex in unique(pDT_ex$celltype)){
    pDT_ex1 <- pDT_ex[celltype == typex]
    for(gx in unique(pDT_ex1[mixscape_class.global != "NTC"]$gene)){
      x1 <- pDT_ex1[gene == gx]$traj
      x2 <- pDT_ex1[gene == "NTC"]$traj
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
  res_ex[, padj.wx := p.adjust(p.wx, method="BH")]
  res_ex[, padj.ks := p.adjust(p.ks, method="BH")]
  write.tsv(res_ex, outS("Statistics.tsv"))
  
  
  # . load ------------------------------------------------------------------
  #res <- fread(outS("Statistics.tsv"))
  
  # . plot stats ------------------------------------------------------------
  ggplot(res_ex, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
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
  ggplot(res_ex, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
    geom_point() +
    theme_bw(12)+ 
    scale_color_gradient2(low="blue", high="red") +
    xRot()
  ggsave(outS("Statistics.pdf"), w=4,h=10)
  
  # . UMAP ------------------------------------------------------------------
  pDT_ex.UMAP <- pDT_ex[abs(traj.scale) < 3]
  ggplot(pDT_ex.UMAP, aes(x=UMAP1, y=UMAP2)) + 
    stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
    theme_bw(12) +
    scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
  ggsave(outS("UMAP.pdf"), w=5,h=5)
  
  
  # . celltypes -------------------------------------------------------------
  ggplot(pDT_ex, aes(x=UMAP1, y=UMAP2, color=celltype)) + 
    geom_point() +
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP.jpg"), w=5,h=5)
  
  ggplot(pDT_ex, aes(x=UMAP1, y=UMAP2)) + 
    geom_hex() +
    facet_grid(. ~ celltype) + 
    theme_bw(12)
  ggsave(outS("Celltypes_UMAP_hex.pdf"), w=20,h=5)
  
  ggplot(pDT_ex, aes(x=celltype)) + 
    geom_bar() +
    scale_y_log10() +
    theme_bw(12)
  ggsave(outS("Celltypes_Nubmers.pdf"), w=5,h=5)
  
  # . plot distributions ex.vivo----------------------------------------------------
  pDT_ex.distr <- copy(pDT_ex)
  ggplot(pDT_ex.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
  pDT_ex.distr <- pDT_ex.distr[abs(traj.scale) < 3]
  pDT_ex.sum <- pDT_ex.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene", "celltype")]
  pDT_ex.stats <- copy(res_ex)
  pDT_ex.stats[, celltype := type]
  pDT_ex.stats[, type := "not.sig"]
  pDT_ex.stats[padj.ks < 0.1, type := "sig.low"]
  pDT_ex.stats[padj.ks < 0.01, type := "sig.high"]
  pDT_ex.distr <- merge(pDT_ex.distr, pDT_ex.stats[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)
  pDT_ex.distr[gene == "NTC", type := "NTC"]
  ggplot(pDT_ex.distr, aes(y=gene, x=traj)) + 
    geom_violin(color=NA, aes(fill=type), scale="width") + 
    scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
    geom_point(data=pDT_ex.sum, color="black") + 
    geom_errorbarh(data=pDT_ex.sum, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
    theme_bw(12) + 
    facet_grid(. ~ celltype, scale="free_x") +
    xRot()
  ggsave(outS("Distribution_ex.vivo.pdf"), w=20,h=10)
  
  ggplot(pDT_ex.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "Mye"],
         aes(x=traj, color=gene)) + 
    theme_bw(12) +
    geom_density()
  ggsave(outS("Distribution_Brd9_ex.vivo.pdf"), w=5,h=4)
  
  
  # . scatterplot -----------------------------------------------------------
  pDT_ex2 <- pDT_ex[celltype %in% c("Ery", "Mye")][, median(traj), by=c("celltype", "gene")]
  pDT_ex2[, V1 := scale(V1), by="celltype"]
  pDT_ex2 <- dcast.data.table(pDT_ex2, gene ~ celltype, value.var = "V1")
  ggplot(pDT_ex2, aes(x=Ery, y=Mye)) + 
    geom_point() + 
    theme_bw(12) +
    geom_abline() + 
    geom_text_repel(aes(label=gene))
  ggsave(outS("Scatter_EryVsMye_ex.vivo.pdf"),w=8,h=8)
}
#test ex.vivo-------------------------------------------
# . test ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res_ex <- data.table()
for(typex in unique(pDT_ex$celltype)){
  pDT1 <- pDT_ex[celltype == typex]
  for(gx in unique(pDT1[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT1[gene == gx]$traj
    x2 <- pDT1[gene == "NTC"]$traj
    if(length(x1) > 10 & length(x2) > 10){
      res_ex <- rbind(res_ex, data.table(
        p.wx=wilcox.test(x1, x2)$p.value,
        p.ks=ks.test(x1, x2)$p.value,
        d=median(x1) - median(x2),
        type=typex,
        gene=gx
      ))
    }
  }
}
res_ex[, padj.wx := p.adjust(p.wx, method="BH")]
res_ex[, padj.ks := p.adjust(p.ks, method="BH")]
write.tsv(res_ex, out("Statistics_ex.vivo.tsv"))
res_ex <- fread(out("Statistics_ex.vivo.tsv"))
ggplot(res_ex, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison_ex.vivo.pdf"),w=20,h=6)

#
head(res_ex)
ggplot(res_ex, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_ex.vivo.pdf"), w=4,h=10)

# . UMAP ------------------------------------------------------------------
pDT.UMAP_ex <- pDT_ex[abs(traj.scale) < 3]
ggplot(pDT.UMAP_ex, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_ex.vivo.pdf"), w=5,h=5)
#

# . plot distributions ----------------------------------------------------
pDT.distr_ex <- copy(pDT_ex)
ggplot(pDT.distr_ex, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~celltype, scales = "free")
pDT.distr_ex <- pDT.distr_ex[abs(traj.scale) < 3]
pDT.sum_ex <- pDT.distr_ex[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("tissue","gene", "celltype")]

pDT.sum_ex[, traj_diff := traj]

cellt<-"Ery"
gen="j"
for (cellt in unique(pDT.sum_ex$celltype)){
  for (gen in unique(pDT.sum_ex$gene)){
    pDT.sum_ex[celltype==cellt & gene == gen]$traj_diff = pDT.sum_ex[celltype==cellt & gene == gen]$traj_diff-pDT.sum_ex[celltype == cellt & gene == "NTC"]$traj}
}

pDT.stats_ex <- copy(res_ex)
pDT.stats_ex[, celltype := type]
pDT.stats_ex[, type := "not.sig"]
pDT.stats_ex[padj.ks <= 0.1, type := "sig.low"]
pDT.stats_ex[padj.ks <= 0.01, type := "sig.high"]
head(pDT.stats_ex)
pDT.distr_ex <- merge(pDT.distr_ex, pDT.stats_ex[,c("gene", "celltype", "type"),with=F], by=c("gene", "celltype"), all.x=TRUE)

pDT.distr_ex[gene=="NTC", type := "NTC"]


ggplot(pDT.distr_ex[gene=="NTC"&celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()

ggplot(pDT.distr_ex, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_ex, color="black") + 
  geom_errorbarh(data=pDT.sum_ex, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype, scale="free_x") +
  xRot()
ggsave(out("Distribution_ex.vivo.pdf"), w=20,h=10)
########################################################
ggplot(pDT.distr_in, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum[gene %in% pDT.distr_ex$gene & celltype =="Mye"], color="black") + 
  geom_errorbarh(data=pDT.sum[gene %in% pDT.distr_ex$gene & celltype =="Mye"], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + xlim(0,40)+
  #facet_grid(. ~ celltype) +
  xRot()
ggsave(outS("Distribution_invivo_Mye.pdf"), w=20,h=10)


ggplot(pDT.distr_ex[celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_ex[gene %in% pDT.distr_ex$gene & celltype =="Mye"], color="black") + 
  geom_errorbarh(data=pDT.sum_ex[gene %in% pDT.distr_ex$gene & celltype =="Mye"], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) +  xlim(0,40)+
  facet_grid(. ~ celltype) +
  xRot()
ggsave(outS("Distribution_exvivo_Mye.pdf"), w=20,h=10)


############################################
#combine-------------------------

pDT.sum_ex[gene=="NTC"]
pDT.sum_combined<-rbind(pDT.sum[gene %in% pDT.distr_ex$gene],pDT.sum_ex)
head(pDT.sum_combined)
pDT.sum_combined_Mye<-pDT.sum_combined[celltype=="Mye"]


Mye<-pDT.distr_ex[celltype=="Mye"][,c("gene","rn","UMAP1","UMAP2","traj","sample.x","celltype","tissue","timepoint","type","traj.scale")]
colnames(Mye)<-gsub("sample.x","sample",colnames(Mye))

pDT.distr_in<-pDT.distr[,c("gene","rn","UMAP1","UMAP2","traj","sample","celltype","tissue","timepoint","type","traj.scale")]
pDT.distr_in<-pDT.distr_in[gene %in% pDT.distr_ex$gene & celltype =="Mye"]
Mye<-rbind(Mye,pDT.distr_in)


#combined_plts-------------------------------

ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_combined_Mye[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.sum_combined_Mye[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_ex.in.pdf"), w=20,h=10)

head(pDT.sum_combined_Mye[gene %in% pDT.distr_ex$gene])
ggplot(data=pDT.sum_combined_Mye[gene %in% pDT.distr_ex$gene], aes(y=gene, x=traj_diff)) + 
  geom_point(data=pDT.sum_combined_Mye[gene %in% pDT.distr_ex$gene], color="black") + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_ex.in_traj-ntc.pdf"), w=20,h=10)

Mye<-Mye[celltype=="Mye"]

#scatter not scaled

pDT2_com <- Mye[tissue %in% c("ex.vivo", "in.vivo")][, median(traj), by=c("tissue", "gene")]

#pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")

ggplot(pDT2_com, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_ex_in.vivo.pdf"),w=8,h=8)

#scatter scaled
pDT2_com <- Mye[tissue %in% c("ex.vivo", "in.vivo")][, median(traj), by=c("tissue", "gene")]
pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")
head(pDT2_com)
ggplot(pDT2_com, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point() + 
  geom_abline()+
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_ex_in.vivo_scaled.pdf"),w=8,h=8)



#traj_median(gene)-ntc

pDT_summary_Mye<-pDT.sum_combined_Mye[,c("tissue","gene","traj_diff")]
head(pDT_summary_Mye)
pDT_summary_Mye<- dcast.data.table(pDT_summary_Mye, gene ~ tissue, value.var = "traj_diff")
pDT_summary_Mye[,diff:=ex.vivo-in.vivo]
head(pDT_summary_Mye)
ggplot(pDT_summary_Mye, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_ex_in.vivo_traj-ntc.pdf"),w=8,h=8)



# logfc correlation and trajectories
col<-merge(pDT_summary_Mye,correlation_logfc_in_ex,by="gene")
col[,correlation_logfc_in_ex:= correlation]
p<-ggplot(col, aes(x=ex.vivo, y=in.vivo))+ 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))

pdf(out("Scatter_ex_in.vivo_traj-ntc_color.pdf"))
p + geom_point(aes(colour = col$correlation))
dev.off()

ggplot(col, aes(x=ex.vivo, y=in.vivo)) + 
  scale_color_manual(col$correlation)+
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_ex_in.vivo_traj-ntc___test.pdf"),w=8,h=8)


#scale to sd(NTC)

Mye[, traj_ntc_norm := traj]
tis="ex.vivo"
gen="j"
for (tis in unique(Mye$tissue)){
  Mye[tissue==tis]$traj_ntc_norm = (Mye[tissue==tis]$traj-(mean(Mye[tissue==tis & gene == "NTC"]$traj)))/sd(Mye[tissue==tis & gene == "NTC"]$traj)
}



###########################
#summary-median of (gene_traj-mean_ntc)/sd(ntc)
pDT.ntc_sum<-  Mye[, .(median_ntc_norm_traj = median(traj_ntc_norm), q1 = quantile(traj_ntc_norm, 0.25), q2 = quantile(traj_ntc_norm, 0.75)),by=c("tissue","gene")]
pDT.ntc_sum[, traj_ntc_norm:= median_ntc_norm_traj]


pDT.ntc_sum_mean<-  Mye[, .(mean_ntc_norm_traj = mean(traj_ntc_norm), q1 = mean(traj_ntc_norm)-sd(traj_ntc_norm), q2 = mean(traj_ntc_norm)+sd(traj_ntc_norm)),by=c("tissue","gene")]
pDT.ntc_sum_mean[, traj_ntc_norm:= mean_ntc_norm_traj]




#median
ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_sum[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_sum[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(out("Distribution_Median_Mye_ex_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)


ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_sum[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_sum[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()

head(pDT.ntc_sum)
pDT.ntc_sum_med<-pDT.ntc_sum[,c("tissue","gene","median_ntc_norm_traj")]
colnames(pDT.ntc_sum_med)<-c("tissue","gene","traj_ntc_norm")
head(pDT.ntc_sum_med)

pDT.ntc_sum_med_diff<- dcast.data.table(pDT.ntc_sum_med, gene ~ tissue, value.var = "traj_ntc_norm")
pDT.ntc_sum_med_diff<-pDT.ntc_sum_med_diff[gene!="Chaf1a"]
pDT.ntc_sum_med_diff[, diff:= in.vivo-ex.vivo]
head(pDT.ntc_sum_med_diff)

ggplot() + 
  #geom_segment(data=pDT.ntc_sum_med_diff, mapping=aes(x=ex.vivo, y=in.vivo, xend=ex.vivo+diff, yend=in.vivo), arrow = arrow(length = unit(0.05, "cm")), size=0.05, color="blue") + 
  geom_point(data=pDT.ntc_sum_med_diff, mapping=aes(x=ex.vivo, y=in.vivo), size=1, shape=21, fill="white") +
  theme_bw(12) +geom_text_repel(aes(label=gene))
ggsave(out("segment.pdf"),w=8,h=8)
head(pDT.ntc_sum_med_diff)
pDT.ntc_sum_med_diff<-reorder(pDT.ntc_sum_med_diff$gene,t)
?reorder
pDT.ntc_sum_med_diff<-pDT.ntc_sum_med_diff[match(t, pDT.ntc_sum_med_diff$gene),]
t
ggplot(pDT.ntc_sum_med_diff, aes(x=ex.vivo, y=gene)) + 
  geom_point(size=1, shape=21, fill="white")+
  geom_segment(data=pDT.ntc_sum_med_diff, mapping=aes(x=ex.vivo, y=gene, 
    xend=ex.vivo+diff, yend=gene),arrow = arrow(length = unit(0.09, "cm")), 
    size=0.09, color="black")+
  theme_bw(12)+
  geom_text_repel(aes(label=gene,max.overlaps =40))

ggsave(out("thenga.png"))
#mean
ggplot(Mye[celltype=="Mye"], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_sum_mean[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_sum_mean[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(out("Distribution_Mean_Mye_ex_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)

##
pDT.ntc_sum_mean<-pDT.ntc_sum_mean[,c("tissue","gene","traj_ntc_norm")]
head(pDT.ntc_sum_mean)

pDT.ntc_sum_mean_diff<- dcast.data.table(pDT.ntc_sum_mean, gene ~ tissue, value.var = "traj_ntc_norm")
pDT.ntc_sum_mean_diff<-pDT.ntc_sum_mean_diff[gene!="Chaf1a"]
pDT.ntc_sum_mean_diff[, diff:= in.vivo-ex.vivo]
head(pDT.ntc_sum_mean_diff)

ggplot(pDT.ntc_sum_med_diff, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point(size=2, shape=21, fill="white") +
  geom_segment(data=pDT.ntc_sum_mean_diff, mapping=aes(x=ex.vivo, 
                                                       y=in.vivo, xend=ex.vivo+diff, yend=in.vivo), 
               arrow = arrow(length = unit(0.2, "cm")), size=0.2, color="black") + 
  theme_bw(12)+ggrepel
ggsave(out("segment.pdf"),w=8,h=8)
########################################
#
pDT.ntc_sum<-pDT.ntc_sum[,c("gene","tissue","median_ntc_norm_traj")]

pDT.ntc_sum_plot<-dcast.data.table(pDT.ntc_sum, gene ~ tissue, value.var = "median_ntc_norm_traj")

pDT.ntc_sum_plot_color<-merge(pDT.ntc_sum_plot,correlation_logfc_in_ex,by="gene")

NTC<-pDT.ntc_sum_plot[gene=="NTC"]
correlation<-0.00
NTC<-cbind(NTC,correlation)

pDT.ntc_sum_plot_color<-rbind(pDT.ntc_sum_plot_color,NTC)
tail(pDT.ntc_sum_plot_color)
#scatter_plot
pDT.ntc_sum_plot[gene=="NTC"]
pDTA<-pDT.ntc_sum_plot_color


#normal-with NTC
ggplot(pDT.ntc_sum_plot, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point()+
  theme_bw(12) +
  geom_text_repel(aes(label=gene))

ggsave(out("Scatter_Mye_ex_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=7,h=8)

#color-No-NTC
ggplot(pDTA, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point(aes(color = pDTA$correlation))+
  
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_ex_in.vivo_traj_ntc_mean_sd_norm_color_corrlogfc.pdf"),w=8,h=8)

#
