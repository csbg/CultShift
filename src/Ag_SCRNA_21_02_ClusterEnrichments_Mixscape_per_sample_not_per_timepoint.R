source("src/00_init.R")
out <- dirout("Ag_SCRNA_21_02_ClusterEnrichments_simple/")


require(doMC)
registerDoMC(10)


# Load annotation ---------------------------------------------------------
inDir <- dirout("Ag_ScRNA_17_proj_celltype_in_monocle_obj")
annList <- fread(inDir("exvivo_Annotations.tsv"))
annList <- rbind(fread(inDir("invivo_Annotations.tsv")),annList)

# Cell types --------------------------------------------------------------
celltype.ann <- CLEAN.CELLTYPES

annList <- annList[!is.na(mixscape_class)]%>%
  filter(mixscape_class.global %in% c("KO","NTC"))%>%
  filter(celltype_projection %in% c(
    "Gran. P", "GMP", "Gran.", "HSC", "MkP", "Mono", "MEP (early)",
    "Eo/Ba")
  )




# Identify and remove or label perturbed clusters (clusters with only perturbations) ---------------------
clDT <- dcast.data.table(annList[,.N, 
                                 by=c("mixscape_class.global",
                                      "tissue", 
                                      "celltype_projection"
                                      )], 
                           tissue + 
                           celltype_projection  ~ mixscape_class.global, value.var = "N")
clDT[, frac := KO/NTC]

# remove clusters from in vivo
(clDT.remove <- clDT[tissue %in% c("ex.vivo","in.vivo"),][(NTC < 5 | is.na(NTC)) & (frac > 25 | is.na(frac))])
clDT
write.tsv(clDT.remove, out("ClustersRemoved.tsv"))
clDT[, remove := ifelse((NTC < 5 | is.na(NTC)) & (frac > 25 | is.na(frac)), "remove", "keep")]

# View the updated table
head(clDT)
clDT <- clDT[remove == "keep"]

annList[tissue == "ex.vivo" & celltype_projection %in% c("Ery","MEP"), Clusters := "remove"]
annList[tissue == "in.vivo" & celltype_projection %in% c("GMP (late)"), Clusters := "remove"]
write.tsv(clDT.remove, out("ClustersRemoved_in.vivo.tsv"))

# Define sets to analyze across tissues ------------------------------------------------------
fish.test.sets <- list()
TISSUES <- unique(annList$tissue)

CELLTYPES <- unique(annList$celltype_projection) 


annList[tissue == tx]
# broad clusters
for(tx in TISSUES){
  x <- annList[tissue == tx]
  x <- merge(x, celltype.ann, by.x="celltype_projection", by.y="Name")[,-"celltype_projection",with=F]
  x[,Clusters := Type]
  fish.test.sets[[paste("broad", tx, sep="_")]] <- x
}
head(fish.test.sets$detailed_ex.vivo)


# Numeric clusters - uses all clusters
for(tx in TISSUES){
  x <- annList[tissue == tx]
  x[, Clusters := celltype_projection]
  fish.test.sets[[paste("detailed", tx, sep="_")]] <- x
}




# Remove (don't use) Mixscape ------------------------------------------------------------
fish.test.sets <- c(
  setNames(lapply(fish.test.sets, 
                  function(dt) dt[mixscape_class.global != "NP"]), 
           paste0(names(fish.test.sets), "_withMixscape")),
  setNames(fish.test.sets, paste0(names(fish.test.sets), "_noMixscape"))
  )



# Remove lowly represented guides -----------------------------------------
sapply(fish.test.sets, nrow)
fish.test.sets <- fish.test.sets[sapply(fish.test.sets, nrow) > 1]



# Perform fisher exact terst ----------------------------------------------
(fish.test.x <- names(fish.test.sets)[1])
foreach(fish.test.x = names(fish.test.sets)) %dopar% {
  res <- data.table()
  pDT1 <- fish.test.sets[[fish.test.x]]
  write.tsv(pDT1[,"rn",with=F], out("Guides_Fisher_Mixscape_",fish.test.x,"_Cells.tsv"))
  stopifnot(sum(is.na(pDT1$CRISPR_Cellranger)) == 0)
  unique(pDT1$sample_broad)
  sx <- pDT1$sample_broad[1]
  for(sx in unique(pDT1$sample_broad)){
    pDT2 <- pDT1[sample_broad == sx]
    gx <- pDT2$genotype[1]
    for(gx in unique(pDT2[mixscape_class != "NTC"]$genotype)){
      cx <- pDT2$Clusters[1]
      for(cx in unique(pDT2$Clusters)){
        ntc <- unique(pDT2[mixscape_class == "NTC"][,.N, by="CRISPR_Cellranger"][N > 20]$CRISPR_Cellranger)[]
        for(ntc in unique(pDT2[mixscape_class == "NTC"][,.N, by="CRISPR_Cellranger"][N > 20]$CRISPR_Cellranger)){
          mx <- as.matrix(with(pDT2[genotype == gx | CRISPR_Cellranger == ntc], table(Clusters == cx, genotype == gx)))
          if(all(dim(mx) == c(2,2))){
            fish <- fisher.test(mx)
            res <- rbind(res, data.table(
              Clusters=cx, 
              mixscape_class=gx, 
              ntc=ntc, 
              p=fish$p.value, 
              OR=fish$estimate, 
              sample=sx, 
              total.cells=sum(mx),
              guide.cells=nrow(pDT2[genotype  == gx])
            ))
          }
        }
      }
    }
  }
  
  # save file
  if(nrow(res) < 3) {
  res[,padj := p.adjust(p, method="BH")]
  res[, log2OR := pmax(-5, pmin(5, log2(OR)))]
  res[,grp := paste(mixscape_class, sample)]
  write.tsv(res[,-c("grp"), with=F], out("Guides_Fisher_Mixscape_",fish.test.x,"per_sample.tsv"))
  
  # plot
  res <- hierarch.ordering(res, toOrder = "grp", orderBy = "Clusters", value.var = "log2OR", aggregate = TRUE)
  ggplot(res, aes(
    x=Clusters,
    y=ntc,
    color=log2OR,
    size=pmin(-log10(padj), 5))) +
    geom_point(shape=16) +
    scale_color_gradient2(name="log2OR", low="blue", high="red") +
    scale_size_continuous(name="padj") +
    facet_grid(mixscape_class + sample ~ ., space = "free", scales = "free") +
    theme_bw(12) +
    theme(strip.text.y = element_text(angle=0)) +
    xRot()
  ggsave(out("Guides_Fisher_Mixscape_",fish.test.x,".pdf"), w=10, h=length(unique(res$grp)) * 0.50 + 1, limitsize = FALSE)
}
}

