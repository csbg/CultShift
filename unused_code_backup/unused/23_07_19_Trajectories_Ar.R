source("src/00_init.R")


library(pheatmap)
require(ggrepel)
require(tidyverse)
ann_in<-fread(dirout_load(paste0("SCRNA_20_Summary/","in.vivo","_monocle.singleR"))("Annotation.tsv"))
ann_ex<-fread(dirout_load(paste0("SCRNA_20_Summary/","ex.vivo","_monocle.singleR"))("Annotation.tsv"))
for (tissue in c("in.vivo","ex.vivo")){
  basedir <- dirout(paste0("SCRNA_50_01_Trajectories/",tissue))
  ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue,"_monocle.singleR"))("Annotation.tsv"))

  # Load annotation ---------------------------------------------------------
  (load(PATHS$SCRNA$MONOCLE.DIR(tissue)))
# ann <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_celltypes.RDS"))
# ann <- ann[rn %in% colnames(monocle.obj)]
  celltypes <- fread('metadata/FIGS_celltypes.tsv')
  celltypes <- celltypes[Name %in% ann$Clusters]

  celltypes_all<-unique(celltypes[Type != "HSC"]$Type)
  celltypes_all
celltypes_all<-c("Mye", "Mega", "Ery", "B-cell", "unclear")
#cds <- monocle.obj[,sample(colnames(monocle.obj), 500)]

#typex <- "Mye"
  for(typex in celltypes_all){


    print(typex)
    cells <- ann[Clusters %in% celltypes[Type %in% c(typex, "HSC")]$Name]$rn
    cells
  #cells <- sample(cells, 100)
    cds <- monocle.obj[,cells]
  
  # recluster (not used but necessary)
    cds <- monocle3::cluster_cells(cds)
  
  # learn graph
   cds <- learn_graph(
      cds,
      verbose = TRUE,
      use_partition = FALSE,
      close_loop = TRUE)
  
  # order cells
  #functional.clusters changed to clusters.final #Ar
    cds <- order_cells(cds, root_cells=ann[rn %in% colnames(cds)][clusters.final == "HSC"]$rn)
  
  
  # Make plot
    tissue
    plot_cells(cds,color_cells_by = "pseudotime")
    
    ggsave(basedir(paste0("Pseudotime_"),typex,".jpg"), w=10,h=10)
  
  # export table
    write.tsv(data.table(data.frame(pseudotime(cds)), keep.rownames = TRUE), out(paste0("Pseudotime_",tissue),typex,".tsv"))
    }
}


# test statistcs ---------------------------------------------------------------
ff <- list.files(out(""), pattern="Pseudotime.*.tsv")
names(ff) <- gsub("Pseudotime_all_(.+).tsv", "\\1", ff)


pDT <- rbindlist(lapply(ff, function(fx) fread(out(fx))), idcol = "celltype")
pDT[, traj := pseudotime.cds.]
pDT$pseudotime.cds. <- NULL

write.tsv(pDT, out("Values_ex_in.tsv"))
#incorporating exvivo and invivo
pDT_in<- merge(pDT, ann_in[Clusters != "HSC"], by="rn")[!is.na(mixscape_class.global)]
pDT_in[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT_in[, traj.scale := scale(traj), by="celltype"]
#only 14d
pDT_in <- pDT_in[timepoint != "28d"]
#
pDT_ex <- merge(pDT, ann_ex[Clusters != "HSC"], by="rn")[!is.na(mixscape_class.global)]
pDT_ex[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT_ex[, traj.scale := scale(traj), by="celltype"]
pDT_ex <- pDT_ex[timepoint != "7d"]
pDT<-rbind(pDT_in,pDT_ex)

common_gene<-unique(pDT_ex$gene)
pDT_common<- pDT[pDT$gene %in% common_gene]
# . test ------------------------------------------------------------------
typex <- "Ery"
gx <- "Rcor1"
res <- data.table()
for(typex in unique(pDT$celltype)){
  pDT3 <- pDT[celltype == typex]
  for(gx in unique(pDT[mixscape_class.global != "NTC"]$gene)){
    x1 <- pDT3[gene == gx]$traj
    x2 <- pDT3[gene == "NTC"]$traj
    if(length(x1) > 10 & length(x2) > 10){
      res <- rbind(res, data.table(
        p.wx=wilcox.test(x1, x2)$p.value,
        p.ks=ks.test(x1, x2)$p.value,
        d=median(x1) - median(x2),
        dmean= mean(x1) -mean(x2),
        
        type=typex,
        gene=gx
      ))
    }
  }
}

################################################
#Ks test to compare the mean(t test in this case because y is also numeric) of the celltypes(eg:Mye) trajectories of a particular KO to NTC
#Wilcox.test ranked sum test between NTC and selected gene KO
res[, padj.wx := p.adjust(p.wx, method="BH")]
res[, padj.ks := p.adjust(p.ks, method="BH")]
write.tsv(res, out("Statistic_ex_in.tsv"))


# . load ------------------------------------------------------------------
res <- fread(out("Statistic_ex_in.tsv"))

res_common<-res[gene %in% common_gene ]

# . plot stats ------------------------------------------------------------
ggplot(res, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison_ex_in.pdf"),w=20,h=6)




# xDT <- melt(res, id.vars = c("type", "gene"))
# xDT[, measurement := gsub("\\..+$", "", variable)]
# xDT[, type := gsub("^.+?\\.", "", variable)]
ggplot(res, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_ex_in.pdf"), w=4,h=10)



# . UMAP ------------------------------------------------------------------
pDT_ex.UMAP <- pDT_ex[abs(traj.scale) < 3]
ggplot(pDT_ex.UMAP, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_ex.pdf"), w=5,h=5)


pDT_in.UMAP <- pDT_in[abs(traj.scale) < 3]
ggplot(pDT_in.UMAP, aes(x=UMAP1, y=UMAP2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_in.pdf"), w=5,h=5)


# . plot distributions ----------------------------------------------------
pDT.distr <- copy(pDT_common)

pDT.distr[, tissue_cell := celltype]
pDT.distr[, tissue := celltype]
pDT.distr[,tissue := gsub("all_","",tissue)]
pDT.distr[,tissue := gsub("in.vivo[a-z0-9A-Z,-]*","in.vivo",tissue)]
pDT.distr[,tissue := gsub("ex.vivo[a-z0-9A-Z,-]*","ex.vivo",tissue)]
pDT.distr[,celltype := gsub("all_in.vivo","",celltype)]
pDT.distr[,celltype := gsub("all_ex.vivo","",celltype)]
ggplot(pDT.distr, aes(x=traj, y=traj.scale)) + geom_hex() + facet_wrap(~tissue_cell, scales = "free")

pDT.distr <- pDT.distr[abs(traj.scale) < 3]
pDT.sum <- pDT.distr[, .(traj = median(traj), q1 = quantile(traj, 0.25), q2 = quantile(traj, 0.75)),by=c("gene","celltype","tissue")]
head(pDT.sum)



pDT.stats[, tissue_cell := type]
pDT.stats[, tissue := type]
pDT.stats[,tissue := gsub("all_","",tissue)]
pDT.stats[,tissue := gsub("in.vivo[a-z0-9A-Z,-]*","in.vivo",tissue)]
pDT.stats[,tissue := gsub("ex.vivo[a-z0-9A-Z,-]*","ex.vivo",tissue)]
pDT.stats[, celltype := tissue_cell]
pDT.stats[,celltype := gsub("all_in.vivo","",celltype)]
pDT.stats[,celltype := gsub("all_ex.vivo","",celltype)]


pDT.stats[, type := "not.sig"]
pDT.stats[padj.ks < 0.1, type := "sig.low"]
pDT.stats[padj.ks < 0.01, type := "sig.high"]
head(pDT.stats)
head(pDT.sum)
head(pDT.distr)
pDT.distr <- merge(pDT.distr, pDT.stats[,c("gene", "tissue_cell", "type"),with=F], by=c("gene", "tissue_cell"), all.x=TRUE)
pDT.distr[is.na(type), type := "NTC"]
head(pDT.distr)

pDT_Mye<-pDT.distr[tissue_cell %in% c("ex.vivoMye","in.vivoMye")]
pDT.sum_Mye<-pDT.sum[celltype %in% c("ex.vivoMye","in.vivoMye")]
ggplot(pDT_Mye, aes(y=gene, x=traj)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.sum_Mye, color="black") +
  geom_errorbarh(data=pDT.sum_Mye, color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ celltype) +
  xRot()
ggsave(out("Distribution_ex_in_Myeloid.pdf"), w=20,h=10)



#


ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "in.vivoMye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggplot(pDT.distr[gene %in% c("Brd9", "Smarcd2", "Smarcd1", 'NTC')][celltype %in% "ex.vivoMye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggsave(out("Distribution_Brd9_ex_in.pdf"), w=5,h=4)


# . scatterplot -----------------------------------------------------------
head(pDT)
pDT5 <- pDT[celltype %in% c("ex.vivoMye", "in.vivoMye")][, median(traj), by=c("celltype", "gene")]

pDT5[, V1 := scale(V1), by="celltype"]
pDT5 <- dcast.data.table(pDT5, gene ~ celltype, value.var = "V1")
ggplot(pDT5, aes(x=ex.vivoMye, y=in.vivoMye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_exVsin_median_traj.pdf"),w=8,h=8)

#d =gene -ntc
head(d_gene_ntc)
d_gene_ntc <- res

d_gene_ntc <- dcast.data.table(d_gene_ntc, gene ~ type, value.var = "d")
d_gene_ntc <- d_gene_ntc[gene %in% pDT_ex$gene]
head(d_gene_ntc)
ggplot(d_gene_ntc, aes(x=ex.vivoMye, y=in.vivoMye)) + 
  geom_point() +  xlim(-10,10)+ylim(-10,10)+
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))

ggsave(out("Scatter_Mye_exVsin_median_traj-ntc.pdf"),w=8,h=8)

cor_celltype<-as.data.frame(d_gene_ntc[,2:7])

rownames(cor_celltype)<-d_gene_ntc$gene

cor_celltype<-cor_celltype[,c(1:4,6)]
cor_celltype<-cor_celltype%>% na.omit()
nrow(cor_celltype)
test<-cor(cor_celltype)

pheatmap(test,cluster_rows = F,cluster_cols = F)

temp_hm_name <- out("pheatmap.png")


temp_hm_name <- paste(deparse(substitute(mat)),".png", sep="")
save_pheatmap_pdf(all_correlation, filename=temp_hm_name)


pDT5[, V1 := scale(V1), by="celltype"]
pDT5 <- dcast.data.table(pDT5, gene ~ celltype, value.var = "V1")
ggplot(pDT5, aes(x=ex.vivoMye, y=in.vivoMye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_exVsin_median_traj.pdf"),w=8,h=8))




pDT5 <- pDT[celltype %in% c("all_ex.vivoEry", "all_in.vivoEry")][, median(traj), by=c("celltype", "gene")]
pDT5[, V1 := scale(V1), by="celltype"]
pDT5 <- dcast.data.table(pDT5, gene ~ celltype, value.var = "V1")
ggplot(pDT5, aes(x=all_ex.vivoEry, y=all_in.vivoEry)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Ery_exVsin_ex.pdf"),w=8,h=8)


pDT5 <- pDT[celltype %in% c("ex.vivoMye", "in.vivoMye")][, median(traj), by=c("celltype", "gene")]
head(pDT)
pDT5[, V1 := scale(V1), by="celltype"]
?dcast.data.table
pDT5 <- dcast.data.table(pDT5, gene ~ celltype, value.var = "V1")
ggplot(pDT5, aes(x=ex.vivoMye, y=in.vivoMye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_exVsin_ex.pdf"),w=8,h=8)




#. res_common statistics--------------------------------------------------
ggplot(res_common, aes(y=-log10(p.wx+1e-10), x=-log10(p.ks+1e-10), color=d)) + 
  geom_point() +
  geom_point(shape=1, color="black") +
  theme_bw(12)+ 
  geom_text_repel(aes(label=paste(gene)), color="black")+
  scale_color_gradient2(low="blue", high="red") +
  facet_grid(. ~ type)
ggsave(out("Statistics_Comparison_ex_in_common_genes.pdf"),w=20,h=6)


ggplot(res_common, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_ex_in_common_genes.pdf"), w=4,h=10)


pDT5 <- pDT[celltype %in% c("ex.vivoMye", "in.vivoMye") & gene %in% pDT_ex$gene][, median(traj), by=c("celltype", "gene")]
head(pDT)
pDT5[, V1 := scale(V1), by="celltype"]
?dcast.data.table
pDT5 <- dcast.data.table(pDT5, gene ~ celltype, value.var = "V1")
ggplot(pDT5, aes(x=ex.vivoMye, y=in.vivoMye)) + 
  geom_point() + 
  theme_bw(12) +
  geom_abline() + 
  geom_text_repel(aes(label=gene))
ggsave(out("Scatter_Mye_exVsin_ex_commom_genes.pdf"),w=8,h=8)

