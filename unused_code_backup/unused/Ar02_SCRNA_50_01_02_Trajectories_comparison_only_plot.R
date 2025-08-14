source("src/00_init.R")
require(umap)
require(ggrepel)
source("src/Ag_Optimized_theme_fig.R")

out1 <- dirout("SCRNA_50_01_Trajectories/")
base.dir <- dirout("SCRNA_50_02_ProjectionTrajectories/")
out <- dirout(base.dir)
##invivo
ann <- fread(dirout_load("SCRNA_20_Summary/in.vivo_monocle.singleR")("Annotation.tsv"))

source("src/FUNC_ProjecTILs_PLUS_Ar1.R")
#exvivo/leukemia

SANN <- fread(PATHS$SCRNA$ANN)
ann1 <- fread(dirout_load("SCRNA_20_Summary/ex.vivo_monocle.singleR")("Annotation.tsv"))
ann2<- fread(dirout_load("SCRNA_20_Summary/leukemia_monocle.singleR")("Annotation.tsv"))

outS<-dirout(paste0(base.dir,"/exvivo_9d_invivo_14d"))
outSLex<-dirout(paste0(base.dir,"/exvivo_9d_leukemia_6d"))
#basedir_corr <- "SCRNA_33_DE_Nebula_testClustering/Output/exvivo_9dinvivo14d/"
#out_corr <- dirout(basedir_corr)
#function

euclidean <- function(a, b) sqrt(sum((a - b)^2))
# test statistcs in.vivo---------------------------------------------------------------
#####################
#numbers per KO
######################################
ann[, gene := gsub("_.+", "", mixscape_class)]
test<- ann%>%group_by(gene)%>%drop_na()%>%summarise(count=n())
n1<-ggplot(test,aes(x=count,y=gene))+geom_col()+FontSize(y.text = 10)+labs(x="leukemia")
ggsave(out("cells_per_guide_invivo.png"),w=12,h=16)

ann1[, gene := gsub("_.+", "", mixscape_class)]
test1<- ann1%>%group_by(gene)%>%drop_na()%>%summarise(count=n())

gene<-c(test$gene[!test$gene %in% test1$gene])
count<-c(rep(0,length(test$gene[!test$gene %in% test1$gene])))
dt<-data.frame(gene,count)
dt<-as_tibble(dt)
test1<-rbind(test1,dt)
n2<-ggplot(test1,aes(x=count,y=gene))+geom_col()+FontSize(y.text = 10)+labs(x="exvivo")
ggsave(out("cells_per_guide_exvivo.png"),w=12,h=16)

ann2[, gene := gsub("_.+", "", mixscape_class)]

test2<- ann2%>%group_by(gene)%>%summarise(count=n())
test2<-test2%>%na.omit()
gene<-c(test$gene[!test$gene %in% test2$gene])
count<-c(rep(0,length(test$gene[!test$gene %in% test2$gene])))
dt<-data.frame(gene,count)
dt<-as_tibble(dt)
test2<-rbind(test2,dt)
n3<-ggplot(test2,aes(x=count,y=gene))+geom_col()+FontSize(y.text = 10)+labs(x="leukemia")
ggsave(out("cells_per_guide_Leukemia.png"),w=12,h=16)
head(test2)
comb<-merge(test,test1,by="gene")
comb<-merge(comb,test2,by="gene")
colnames(comb)<-c("gene","count_invivo","count_exvivo","count_leukemia")
write.tsv(comb,out("counts_per_KO.tsv"))
n1+n2+n3
ggsave(out("cells_per_guide_per_tissue.png"),w=12,h=16)
###########################
######

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
typex <- "Ery"
gx <- "Rcor1"
UN1Mean<-0
UN1Med<-0
UN2Mean<-0
UN2Med<-0
res_umap <- data.table()
for(typex in unique(pDT$celltype)){
  pDT1 <- pDT[celltype == typex]
  for(gx in unique(pDT1$gene)){
    U1 <- pDT1[gene == gx]$UMAP1
    UN1 <- pDT1[gene == "NTC"]$UMAP1
    U2<-pDT1[gene == gx]$UMAP2
    UN2 <- pDT1[gene == "NTC"]$UMAP2
    res_umap <- rbind(res_umap, data.table(
    U1Med=median(U1),
    U2Med=median(U2),
    UN1Med=median(UN1),
    UN2Med=median(UN2),
    #UN1Mean=mean(UN1),
    #UN2Mean=mean(UN2),
    U1_norm=round(median((U1-mean(UN1))/sd(UN1)),3),
    U2_norm=round(median((U2-mean(UN2))/sd(UN2)),3),
    q11 = quantile(U1, 0.25),
    q12 = quantile(U1, 0.75),
    q21 = quantile(U2, 0.25),
    q22 = quantile(U2, 0.75),
    d=sqrt(sum((median(U1)-median(UN1))^2,(median(U2)-median(UN2))^2)),
    
    type=typex,
    gene=gx
      ))
    
  }
}
res_umap[,tissue:="in.vivo"]
head(pDT)



UMAP <- pDT[abs(traj.scale) < 3]
head(pDT.UMAP)
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

ggplot(pDT.distr[gene %in% c("Brd9", 'NTC')][celltype %in% "Mye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
ggsave(out1("Distribution_Brd9_in.vivo.pdf"), w=5,h=4)


# Annotation ex.vivo--------------------------------------------------------------



ff1 <- list.files(out(""), pattern="Output.*.tsv")
ff1 <- ff1[!grepl("in.vivo", ff1)]
ff1 <- ff1[!grepl("leukemia", ff1)]
names(ff1) <- gsub("Output_(.+).tsv", "\\1", ff1)
ff1
pDT_ex <- rbindlist(lapply(ff1, function(fx) fread(out(fx))), idcol = "sample")
pDT_ex[, traj := pseudotime]
pDT_ex$pseudotime <- NULL
head(pDT_ex)
ann_3<-rbind(ann1,ann2,fill=T)
#write.tsv(pDT_ex, out("Values_ex.vivo.tsv"))
pDT_ex <- merge(pDT_ex, ann_3, by="rn")[!is.na(mixscape_class.global)]
pDT_ex[, gene := gsub("_.+$", "", CRISPR_Cellranger)]
pDT_ex[, celltype := ct]
pDT_ex[, traj.scale := scale(traj), by="celltype"]
pDT_ex[, id := paste0(tissue, "_", timepoint)]
pDT_ex<- pDT_ex[timepoint != "7d"]

#test ex.vivo-------------------------------------------
# . test ------------------------------------------------------------------
for(idx in unique(pDT_ex$id)){
  #outS <- dirout(paste0(base.dir, idx))
  pDT_ex <- pDT_ex[id == idx]
  
  # . test ------------------------------------------------------------------
  typex <- "Ery"
  gx <- "Rcor1"
  res_ex <- data.table()
  for(typex in unique(pDT_ex$celltype)){
    pDT_ex1 <- pDT_ex[celltype == typex]
    for(gx in unique(pDT_ex1[mixscape_class.global != "NTC"]$gene)){
      x1 <- pDT_ex1[gene == gx]$traj
      x2 <- pDT_ex1[gene == "NTC"]$traj
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
  }}
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
ggplot(res_ex, aes(y=gene, x=type, size=pmin(5, -log10(padj.wx)), color=d)) + 
  geom_point() +
  theme_bw(12)+ 
  scale_color_gradient2(low="blue", high="red") +
  xRot()
ggsave(out("Statistics_ex.vivo.pdf"), w=4,h=10)

# . UMAP ------------------------------------------------------------------

typex <- "Ery"
gx <- "Rcor1"
UN1Mean<-0
UN1Med<-0
UN2Mean<-0
UN2Med<-0
res_umap_ex <- data.table()
for(typex in unique(pDT_ex$celltype)){
  pDT1_ex <- pDT_ex[celltype == typex]
  for(gx in unique(pDT1_ex$gene)){
    U1 <- pDT1_ex[gene == gx]$UMAP_1
    UN1 <- pDT1_ex[gene == "NTC"]$UMAP_1
    U2<-pDT1_ex[gene == gx]$UMAP_2
    UN2 <- pDT1_ex[gene == "NTC"]$UMAP_2
    res_umap_ex <- rbind(res_umap_ex, data.table(
      U1Med=median(U1),
      U2Med=median(U2),
      UN1Med=median(UN1),
      UN2Med=median(UN2),
      U1_norm=round(median((U1-mean(UN1))/sd(UN1)),3),
      U2_norm=round(median((U2-mean(UN2))/sd(UN2)),3),
      q11 = quantile(U1, 0.25),
      q12 = quantile(U1, 0.75),
      q21 = quantile(U2, 0.25),
      q22 = quantile(U2, 0.75),
      d=sqrt(sum((median(U1)-median(UN1))^2,(median(U2)-median(UN2))^2)),
      
      type=typex,
      gene=gx
    ))
    
  }
}
res_umap_ex[,tissue:="ex.vivo"]
pDT.UMAP_ex <- pDT_ex[abs(traj.scale) < 3]
ggplot(pDT.UMAP_ex, aes(x=UMAP_1, y=UMAP_2)) + 
  stat_summary_hex(bins = 100, aes(z=traj.scale),fun=mean) +
  theme_bw(12) +
  scale_fill_gradientn(colors=c("lightgrey", "blue", "purple", "red", "orange")) 
ggsave(out("UMAP_ex.vivo.pdf"), w=5,h=5)
#

# . plot distributions ex.vivo----------------------------------------------------
pDT.distr_ex <- copy(pDT_ex)

ggplot(pDT.distr_ex, aes(x=traj, y=traj.scale)) + geom_hex() +
  facet_wrap(~celltype, scales = "free")
pDT.distr_ex <- pDT.distr_ex[abs(traj.scale) < 3]
pDT.sum_ex <- pDT.distr_ex[, .(traj = median(traj), 
                               q1 = quantile(traj, 0.25),
                               q2 = quantile(traj, 0.75)),by=c("tissue","gene", "celltype")]

pDT.sum_ex[, traj_diff := traj]
head(pDT.sum_ex)

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

#################################################################
ggplot(pDT.distr_ex[gene %in% c("Brd9", 'NTC')][celltype %in% "Mye"],
       aes(x=traj, color=gene)) + 
  theme_bw(12) +
  geom_density()
pDT.distr_ex[gene %in% c("Brd9", 'NTC')][celltype %in% "Mye"]
pDT.distr[gene %in% c("Brd9", 'NTC')][celltype %in% "Mye"]
################
############################################
#combine ex.vivo and in.vivo-------------------------

pDT.sum_ex[gene=="NTC"]
pDT.sum_combined<-rbind(pDT.sum[gene %in% pDT.distr_ex$gene],pDT.sum_ex)
pDT.sum_combined_Mye<-pDT.sum_combined[celltype=="Mye"]
pDT.sum_combined_Mye %>% write_rds(base.dir("pDT.sum_combined_Mye.rds"))
Mye<-pDT.distr_ex[celltype=="Mye"][,c("gene","rn","UMAP1","UMAP2",
                                      "traj","sample.x","celltype",
                                      "tissue","timepoint","type","traj.scale")]
colnames(Mye)<-gsub("sample.x","sample",colnames(Mye))

pDT.distr_in <- pDT.distr[,c("gene","rn","UMAP1","UMAP2","traj","sample","celltype","tissue","timepoint","type","traj.scale")]
#only selecting genes also present in exvivo
#pDT.distr_in<-pDT.distr_in[gene %in% pDT.distr_ex$gene & celltype =="Mye"]
Mye <- rbind(Mye,pDT.distr_in)

Mye <- rbind(Mye,pDT.distr_in) 
Mye %>%  write_rds(base.dir("Myeloid_traj.distr.ex.in.rds"))
#combined_plts-------------------------------





Mye<-Mye[celltype=="Mye"]

#scatter not scaled

pDT2_com <- Mye[tissue %in% c("ex.vivo", "in.vivo")][, median(traj), by=c("tissue", "gene")]

#pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")

ggplot(pDT2_com, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point() + 
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(outS("Scatter_ex_in.vivo.pdf"),w=8,h=8)

#scatter scaled
pDT2_com <- Mye[tissue %in% c("ex.vivo", "in.vivo")][, median(traj), by=c("tissue", "gene")]
pDT2_com[, V1 := scale(V1), by="tissue"]
pDT2_com <- dcast.data.table(pDT2_com, gene ~ tissue, value.var = "V1")


#traj_median(gene)-ntc

pDT_summary_Mye<-pDT.sum_combined_Mye[,c("tissue","gene","traj_diff")]

pDT_summary_Mye<- dcast.data.table(pDT_summary_Mye, gene ~ tissue, value.var = "traj_diff")
pDT_summary_Mye[,diff:=ex.vivo-in.vivo]


#scale to sd(NTC)

Mye[, traj_ntc_norm := traj]
tis="ex.vivo"
gen="j"
for (tis in unique(Mye$tissue)){
  Mye[tissue==tis]$traj_ntc_norm = (Mye[tissue==tis]$traj-(mean(Mye[tissue==tis & gene == "NTC"]$traj)))/sd(Mye[tissue==tis & gene == "NTC"]$traj)
}



###########################
#summary-median of (gene_traj-mean_ntc)/sd(ntc)
pDT.ntc_sum_ex<-  Mye[, .(median_ntc_norm_traj = median(traj_ntc_norm), q1 = quantile(traj_ntc_norm, 0.25), q2 = quantile(traj_ntc_norm, 0.75)),by=c("tissue","gene")]
pDT.ntc_sum_ex[, traj_ntc_norm:= median_ntc_norm_traj]


pDT.ntc_sum_ex_mean<-  Mye[, .(mean_ntc_norm_traj = mean(traj_ntc_norm),
                               q1 = mean(traj_ntc_norm)-sd(traj_ntc_norm), q2 = mean(traj_ntc_norm)+sd(traj_ntc_norm)),by=c("tissue","gene")]
pDT.ntc_sum_ex_mean[, traj_ntc_norm:= mean_ntc_norm_traj]




#median
ggplot(Mye[celltype=="Mye" & gene %in% pDT.distr_ex$gene], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_sum_ex[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_sum_ex[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_Median_Mye_ex_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)




pDT.ntc_sum_ex_med<-pDT.ntc_sum_ex[,c("tissue","gene","median_ntc_norm_traj")]
colnames(pDT.ntc_sum_ex_med)<-c("tissue","gene","traj_ntc_norm")


pDT.ntc_sum_ex_med_diff<- dcast.data.table(pDT.ntc_sum_ex_med, gene ~ tissue,
                                           value.var = "traj_ntc_norm")

pDT.ntc_sum_ex_med_diff[, diff:= in.vivo-ex.vivo]

pDT.ntc_sum_ex_med_diff <- pDT.ntc_sum_ex_med_diff %>%na.omit()
pDT.ntc_sum_ex_med_diff %>%
  write_rds(base.dir("pDT.ntc_sum_ex_med_diff.rds"))
#################
#paper fig --------------------
#################
means <- Mye %>%
  filter(gene %in% c("Spi1", "NTC")) %>%
  group_by(gene, tissue) %>%
  summarize(mean_traj = mean(traj_ntc_norm, na.rm = TRUE))

# Plot density with mean lines

density <- ggplot(Mye %>% filter(gene %in% c("Spi1", "NTC")),
       aes(x = traj_ntc_norm, color = gene)) + 
  facet_grid(cols = vars(tissue)) +
  theme_bw(12) +
  geom_density() +
  geom_vline(data = means, aes(xintercept = mean_traj, color = gene), linetype = "dashed") +  # Add mean lines
  optimized_theme_fig() +
  theme(legend.position = "right")+labs(x=NULL)
ggsave(outS("density.pdf"))

distr <- ggplot(pDT.ntc_sum_ex_med_diff, aes(x=ex.vivo, y=gene)) + 
  geom_point(size=1, shape=21, fill="white")+
  geom_segment(data=pDT.ntc_sum_ex_med_diff, mapping=aes(x=ex.vivo, y=gene, 
                                                      xend=ex.vivo+diff, yend=gene),arrow = arrow(length = unit(0.09, "cm")), 
               size=0.2, color="black")+
  labs(x="(ex-vivo to in-vivo)", y="KOs")+
  theme_bw(12)+optimized_theme_fig()+coord_flip()

ggsave(outS("exvivo_invivo_segment_norm_ntc.pdf"), w=5, h=length(unique(pDT.ntc_sum_ex_med_diff $gene))*.21)
library(patchwork)
combined <- density + distr + plot_layout(heights = c(1,1.5))
basedir<-dirout("Figure2")
ggsave(basedir("fig2.1.pdf"),combined, w=18, h=8, units = "cm")
#mean
ggplot(Mye[celltype=="Mye" & gene %in% pDT.distr_ex$gene], aes(y=gene, x=traj_ntc_norm)) + 
  geom_violin(color=NA, aes(fill=type), scale="width") + 
  scale_fill_manual(values=c(not.sig = "grey", sig.low="#fb9a99", sig.high="#e31a1c", NTC="#a6cee3")) +
  geom_point(data=pDT.ntc_sum_ex_mean[gene %in% pDT.distr_ex$gene], color="black") + 
  geom_errorbarh(data=pDT.ntc_sum_ex_mean[gene %in% pDT.distr_ex$gene], color="black", aes(xmin=q1, xmax=q2), height = .2) + 
  theme_bw(12) + 
  facet_grid(. ~ tissue) +
  xRot()
ggsave(outS("Distribution_Mean_Mye_ex_in.vivo_traj_ntc_mean_sd_norm.pdf"),w=8,h=8)

##
pDT.ntc_sum_ex_mean<-pDT.ntc_sum_ex_mean[,c("tissue","gene","traj_ntc_norm")]
head(pDT.ntc_sum_ex_mean)

pDT.ntc_sum_ex_mean_diff<- dcast.data.table(pDT.ntc_sum_ex_mean, gene ~ tissue, value.var = "traj_ntc_norm")
#pDT.ntc_sum_ex_mean_diff<-pDT.ntc_sum_ex_mean_diff[gene!="Chaf1a"]
pDT.ntc_sum_ex_mean_diff[, diff:= in.vivo-ex.vivo]
head(pDT.ntc_sum_ex_mean_diff)


########################################
#
pDT.ntc_sum_ex<-pDT.ntc_sum_ex[,c("gene","tissue","median_ntc_norm_traj")]

pDT.ntc_sum_ex_plot<-dcast.data.table(pDT.ntc_sum_ex, gene ~ tissue, value.var = "median_ntc_norm_traj")

pDT.ntc_sum_ex_plot_color<-merge(pDT.ntc_sum_ex_plot,correlation_logfc,by="gene")

NTC<-pDT.ntc_sum_ex_plot[gene=="NTC"]
correlation<-0.00
NTC<-cbind(NTC,correlation)

pDT.ntc_sum_ex_plot_color<-rbind(pDT.ntc_sum_ex_plot_color,NTC)
tail(pDT.ntc_sum_ex_plot_color)
#scatter_plot
pDT.ntc_sum_ex_plot[gene=="NTC"]
pDTA<-pDT.ntc_sum_ex_plot_color



#color--NTC
ggplot(pDTA, aes(x=ex.vivo, y=in.vivo)) + 
  geom_point(aes(color = pDTA$correlation))+
  
  theme_bw(12) +
  geom_text_repel(aes(label=gene))
ggsave(outS("Scatter_Mye_ex_in.vivo_traj_ntc_mean_sd_norm_color_corrlogfc.pdf"),w=8,h=8)

####
#UMAP_combined-------------
res_umap_in_ex<-rbind(res_umap,res_umap_ex)
colnames(res_umap_ex)<-paste("ex", colnames(res_umap_ex), sep = "_")
colnames(res_umap_ex)<-gsub("ex_gene","gene",colnames(res_umap_ex))
colnames(res_umap_ex)<-gsub("ex_type","type",colnames(res_umap_ex))
head(res_umap_ex)
res_umap_in_ex_merged<-merge(res_umap,res_umap_ex,by=c("gene","type"))
ggplot(res_umap_in_ex,aes(U1Med,U2Med))+geom_point()+facet_grid(~tissue,scale="free_x")+
  geom_text_repel(aes(label=gene))


ggplot(res_umap_in_ex_merged, aes(x=ex_U1_norm, y=ex_U2_norm)) + 
  geom_point(size=1, shape=21, fill="white")+
  geom_segment(data=res_umap_in_ex_merged, mapping=aes(x=ex_U1_norm, y=ex_U2_norm, 
               xend=ex_U1_norm+(ex_U1_norm-U1_norm), yend=ex_U2_norm+(ex_U2_norm-U2_norm)),
               arrow = arrow(length = unit(0.09, "cm")), 
               size=0.09, color="black")+xlim(-2,2)+geom_text_repel(aes(label=gene),max.overlaps = 40)+
  xlim(-1,2)+ylim(-7,7)+
  theme_bw(12)
ggplot(res_umap_in_ex_merged, aes(x=ex_U1_norm, y=gene)) + 
  geom_point(size=1, shape=21, fill="white")+
  geom_segment(data=res_umap_in_ex_merged, mapping=aes(x=ex_U1_norm, y=gene, 
                                                       xend=ex_U1_norm+(ex_U1_norm-U1_norm), yend=gene),
               arrow = arrow(length = unit(0.09, "cm")), 
               size=0.09, color="black")+facet_grid(~type,scale="free_x")+
  theme_bw(12)+FontSize(y.text = 5)
  ggsave(out_ex_in("UMAP1_norm_NTC_in_ex.png"))
  
  ggplot(res_umap_in_ex_merged, aes(x=ex_U2_norm, y=gene)) + 
    geom_point(size=1, shape=21, fill="white")+
    geom_segment(data=res_umap_in_ex_merged, mapping=aes(x=ex_U2_norm, y=gene, 
                                                         xend=ex_U2_norm+(ex_U2_norm-U2_norm), yend=gene),
                 arrow = arrow(length = unit(0.09, "cm")), 
                 size=0.09, color="black")+facet_grid(~type,scale="free_x")+
    theme_bw(12)+FontSize(y.text = 5)
  ggsave(out_ex_in("UMAP2_norm_NTC_in_ex.png"))
  #euclidean distance
  ggplot(res_umap_in_ex_merged, aes(x=ex_d, y=gene)) + 
    geom_point(size=1, shape=21, fill="white")+
    geom_segment(data=res_umap_in_ex_merged, mapping=aes(x=ex_d, y=gene, 
                                                         xend=d, yend=gene),
                 arrow = arrow(length = unit(0.09, "cm")), 
                 size=0.09, color="black")+facet_grid(~type,scale="free_x")+
    theme_bw(12)+FontSize(y.text = 5)
  ggsave(out_ex_in("euclidean_dist_to_NTC_in_ex.png"))
##############################
