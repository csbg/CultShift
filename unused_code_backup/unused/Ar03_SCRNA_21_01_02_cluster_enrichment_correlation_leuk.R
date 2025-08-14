
source("src/00_init.R")

base.dir <- "SCRNA_50_02_ProjectionTrajectories/"

# . Cluster enrichment analyses ---------------------------------------------
tx <- "ex.vivo"
inDir <- dirout_load("SCRNA_21_02_ClusterEnrichments_simple")
#for(tx in names(inDir.funcs)){
# out directory
out <- dirout(paste0(base.dir, "/", "cluster.enrichments/"))
exvivo_gene<-unique(fread(inDir("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape.tsv"))$mixscape_class)
exvivo_celltype<-unique(fread(inDir("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape.tsv"))$Clusters)
# Collect enrichment scores
typex<-""
gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_basic.*.pdf"))
for(typex in gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_basic.*.pdf"))){
  fish.file <- inDir("Guides_Fisher_Mixscape_basic_leukemia_noMixscape.tsv")
  if(!file.exists(fish.file)) next
  
  fish.full <- fread(fish.file)
  unique(fish.full$mixscape_class)
  
  fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]
  
  
  #fish.full <- merge(fish.full, unique(SANN[,c("sample_broad", "timepoint"),with=F]), by.x="sample", by.y="sample_broad")
  timex <- "9d"
  
  for(timex in c(unique(fish.full$sample))){
    fish <- copy(fish.full)
    typex
    timex
    #fish <- fish[mixscape_class == "Rbbp4"]
    if(timex != "all") fish <- fish[sample == timex]
    
    # summarize across NTCs
    fish <- fish[, .(
      log2OR=mean(log2OR), 
      dir=length(unique(sign(log2OR[padj < 0.01]))) <= 1, 
      #dir=length(unique(sign(log2OR)))==1, 
      padj=sum(padj < 0.01),
      N=.N,
      Ncells=sum(unique(guide.cells))), by=c("sample", "Clusters", "mixscape_class")]
    fish[dir == FALSE, padj := 0]
    fish[dir == FALSE, log2OR := NA]
    
    # legacy
    fish[, gene := mixscape_class]
    
    # Summarize across guides
    # fish[, gene := gsub("_.+", "", mixscape_class)]
    # fish[gene == "Pu.1", gene := "Spi1"]
    # fish <- fish[, .(
    #   log2OR=mean(log2OR, na.rm=TRUE), 
    #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
    #   padj=sum(padj), 
    #   N=sum(N)), by=c("sample", "Clusters", "gene")]
    # fish[dir == FALSE, padj := 0]
    # fish[dir == FALSE, log2OR := NA]
    
    # summarize across samples
    # fish <- fish[, .(
    #   log2OR=mean(log2OR, na.rm=TRUE), 
    #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
    #   padj=sum(padj), 
    #   N=sum(N)), by=c("Clusters", "gene")]
    # fish[dir == FALSE, padj := 0]
    # fish[dir == FALSE, log2OR := NA]
    head(fish)
    #
    fish<-fish[gene %in% exvivo_gene]
    fish<-fish[Clusters %in% exvivo_celltype]
    # setup for plotting
    fish[padj == 0 | is.na(log2OR), log2OR := 0]
    fish[, sig.perc := padj / N]
    fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
    fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "Clusters", value.var = "log2OR")
    fish[, Clusters := gsub("^Gran$", "Gran.", Clusters)]
    #fish[, Clusters := cleanCelltypes(Clusters)]
    #fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")
    ggplot(fish, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
      themeNF(rotate=TRUE) +
      scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
      scale_size_continuous(name="% sign.", range = c(0,5)) +
      geom_point() +
      geom_point(shape=1, color="lightgrey") +
      xlab("Gene") + ylab("Cell type")
    ggsaveNF(
      out("Cluster_enrichments_basic_noMixscape",typex,"_", timex, "Filtered.pdf"), 
      w=length(unique(fish$gene))*0.05 + 0.5,
      h=length(unique(fish$Clusters))*0.05 + 0.5)
    write.tsv(fish, out("Cluster_enrichments_Filtered_basic_noMixscape",typex,"_", timex,".tsv"))
  }
}
#}
# Change to broad or basic celltype classification as required

ex.vivo_9d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_Filtered_basic_noMixscapebasic_ex.vivo_noMixscape_9d.tsv"))
leukemia_6d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_Filtered_basic_noMixscapebasic_leukemia_noMixscape_6d.tsv"))

leukemia_6d<-leukemia_6d %>%filter(Clusters %in% unique(ex.vivo_9d$Clusters)) %>% filter(gene %in% unique(ex.vivo_9d$gene)) %>% filter(Clusters %in% exvivo_celltype)
ex.vivo_9d<-ex.vivo_9d %>%filter(Clusters %in% unique(leukemia_6d$Clusters)) %>% filter(gene %in% unique(leukemia_6d$gene)) %>% filter(Clusters %in% exvivo_celltype)

leukemia_6d_arrange<- leukemia_6d %>% select(gene,Clusters,log2OR) %>%
  group_by(gene, Clusters)
leukemia_6d_arrange<-as.data.table(leukemia_6d_arrange)
leukemia_6d<- leukemia_6d_arrange[with(leukemia_6d_arrange,order(gene,Clusters))]
head(leukemia_6d)

ex.vivo_9d_arrange<- ex.vivo_9d %>% select(gene,Clusters,log2OR) %>% 
  group_by(gene, Clusters)
ex.vivo_9d_arrange<-as.data.table(ex.vivo_9d_arrange)
ex.vivo_9d<- ex.vivo_9d_arrange[with(ex.vivo_9d_arrange,order(gene,Clusters))]


combined_cluster_leuk<-merge(ex.vivo_9d,leukemia_6d,by=c("gene","Clusters")) %>% pivot_wider(
  names_from = gene,
  values_from=c(log2OR.x,log2OR.y)
)

colnames(combined_cluster_leuk)[grep("log2OR.x",colnames(combined_cluster_leuk))]<-gsub("log2OR.x","exvivo",colnames(combined_cluster_leuk)[grep("log2OR.x",colnames(combined_cluster_leuk))])
colnames(combined_cluster_leuk)[grep("log2OR.y",colnames(combined_cluster_leuk))]<-gsub("log2OR.y","leukemia",colnames(combined_cluster_leuk)[grep("log2OR.y",colnames(combined_cluster_leuk))])
unique(ex.vivo_9d$Clusters
       )
unique(leukemia_6d$Clusters)
corr<-"pearson"
#for (corr in c("pearson","spearman")){  
  test2<-cor(combined_cluster_leuk[grep("exvivo",colnames(combined_cluster_leuk))],y=combined_cluster_leuk[grep("leukemia",colnames(combined_cluster_leuk))],
            method = "pearson")
  head(test2)
  pdf(out("Basic_pheatmap_across_all","exvivo9d","exvivo14d",corr,"basic.pdf"))  
  all_correlation<-pheatmap(test2,cluster_rows = F,cluster_cols = F)
  dev.off()
  all_correlation
  
  
  
  diag_correlation_leuk_cluster<-data.frame(test2[row(test2)==col(test2)])
  head(diag_correlation_leuk_cluster)
  colnames(diag_correlation_leuk_cluster)<-"correlation"
  rownames(diag_correlation_leuk_cluster)<-gsub("exvivo_","",rownames(test2))
  ggplot(corr_enrichment_vs_logfc, aes(x=enrichment, y=logfc)) + 
    geom_point() + 
    theme_bw(12) +
    geom_text_repel(aes(label=gene))
  ggsave(out1(paste0("Scatter_enrichment",corr,"_vs_logfc","_nomixscape.pdf")),w=8,h=8)
#}

diag_correlation_leuk_cluster["gene"]<-rownames(diag_correlation_leuk_cluster)
head(diag_correlation_leuk_selected)
colnames(diag_correlation_leuk_cluster)<-c("corr_cluster","gene")
head(ex.vivo_9d)
p<-diag_correlation_leuk_cluster[order(diag_correlation_leuk_cluster$gene), ]
p1<-ggplot(p, aes(x=gene, y=corr_cluster)) +
  geom_bar(stat='identity') + ylim(-1.000,1.000)+labs(y="Cluster-Correlation")+
  theme_bw(12)+FontSize(y.text = 10)+
  coord_flip()+scale_x_discrete(limits = rev)+theme(axis.title.y = element_blank())+FontSize(y.text = 20,x.text = 20)

ggsave(out("correlation_cluster_alphabetic_exvivo_leuk.png"),w=8,h=12)
?coord_flip
