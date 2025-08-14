source("src/00_init.R")
library(tidyselect)
library(pheatmap)
library(cowplot)
library(patchwork)
source("src/00_init.R")
source("Ar02_SCRNA_50_01_02_Trajectories_comparison_only_plot.R")
source("Ar03_SCRNA_21_01_02_cluster_enrichment_correlation_leuk.R")
source("Ar03_SCRNA_21_01_02_cluster_enrichment_correlation.R")

basedir <- "Ar_05_combined_Clustering/Output"
inp<-"SCRNA_33_DE_Nebula_testClustering/"
out <- dirout(basedir)
ex<-"9d"
inv<-"14d"
leuk<-"6d"
#logFC
#for (ex in c("7d","9d")){
 # for (inv in c("14d","28d")){
    exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
        
    invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
        
    leukemia=fread(dirout(paste0(inp, "leukemia_",leuk,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    
    
genes<-unique(exvivo$gene_id)
  combined<-merge(exvivo,invivo,by=c("guide","gene_id")) %>% pivot_wider(
  names_from = guide,
  values_from=c(estimate.x,estimate.y))
  colnames(combined)[grep("estimate.x",colnames(combined))]<-gsub("estimate.x","exvivo",colnames(combined)[grep("estimate.x",colnames(combined))])
  colnames(combined)[grep("estimate.y",colnames(combined))]<-gsub("estimate.y","invivo",colnames(combined)[grep("estimate.y",colnames(combined))])
    
  combined_leuk<-merge(exvivo,leukemia,by=c("guide","gene_id")) %>% pivot_wider(
  names_from = guide,
  values_from=c(estimate.x,estimate.y)
    )
    colnames(combined_leuk)[grep("estimate.x",colnames(combined_leuk))]<-gsub("estimate.x","exvivo",colnames(combined_leuk)[grep("estimate.x",colnames(combined_leuk))])
    colnames(combined_leuk)[grep("estimate.y",colnames(combined_leuk))]<-gsub("estimate.y","leukemia",colnames(combined_leuk)[grep("estimate.y",colnames(combined_leuk))])
    
    
    
    out_ex_in <- dirout(paste0(basedir,"/exvivo_",ex,"invivo",inv))
    out_ex_leuk<- dirout(paste0(basedir,"/exvivo_",ex,"leukemia",leuk))
    corr<-"pearson"
    for (corr in c("spearman","pearson")){
      
      test<-cor(combined[grep("exvivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
                method = corr)
      
      test_leuk<-cor(combined_leuk[grep("exvivo",colnames(combined_leuk))],y=combined_leuk[grep("leukemia",colnames(combined_leuk))],
                method = corr)
      head(test_leuk)
      png(out_ex_in(paste0("all_correlations_","invivo",inv,"exvivo",ex,corr,"_logfc_.png")))
      all_correlation<-pheatmap(test,cluster_rows = T,cluster_cols = F)
      dev.off()
      
      diag_correlation<-data.frame(test[row(test)==col(test)])
      colnames(diag_correlation)<-"correlation"
      rownames(diag_correlation)<-gsub("exvivo_","",rownames(test))
      
      diag_correlation["gene"]<-rownames(diag_correlation)
      
      diag_correlation_selected<-diag_correlation%>% filter(gene %in% unique(exvivo$gene_id))%>% arrange(desc(correlation))
      
      diag_correlation_leuk<-data.frame(test_leuk[row(test_leuk)==col(test_leuk)])
      colnames(diag_correlation_leuk)<-"correlation"
      rownames(diag_correlation_leuk)<-gsub("exvivo_","",rownames(test_leuk))
      
      diag_correlation_leuk["gene"]<-rownames(diag_correlation_leuk)
      head(diag_correlation_leuk)
      diag_correlation_leuk_selected<-diag_correlation_leuk%>% filter(gene %in% unique(exvivo$gene_id))%>% arrange(desc(correlation))
      
      
      write.tsv(data.table(data.frame(diag_correlation), keep.rownames = TRUE),
                out_ex_in(paste0("correlation_invivo",inv,"vs_exvivo",ex,"_logfc_",corr,"_.tsv")))
      write.tsv(data.table(data.frame(diag_correlation_leuk), keep.rownames = TRUE),
                out_ex_in(paste0("correlation_exvivo",inv,"vs_leuk",ex,"_logfc_",corr,"_.tsv")))
      ggplot(diag_correlation, aes(x=reorder(gene,correlation), y=correlation)) +
        geom_bar(stat='identity') + ylim(-0.1000,0.5000)+
        coord_flip()
      
      ggsave(out_ex_in("correlation_exvivo_invivo_logfc_all.png"), w=20,h=10)
      g1<-ggplot(diag_correlation_selected, aes(x=reorder(gene,correlation), y=correlation)) +
        geom_bar(stat='identity') + ylim(-0.1000,0.5000)+FontSize(y.text = 10)+
                coord_flip()+labs( y ="LogFC_correlation", x = "genes")
     
      ggsave(out_ex_in("correlation_exvivo_invivo_logfc_selected.png"), w=20,h=10)
      
      g11<-ggplot(diag_correlation_leuk_selected, aes(x=reorder(gene,correlation), y=correlation)) +
        geom_bar(stat='identity') + ylim(-0.1000,0.5000)+FontSize(y.text = 10)+
        coord_flip()+labs( y ="LogFC_correlation", x = "genes")
      ggsave(out_ex_leuk("correlation_exvivo_leukemia_logfc_selected.png"), w=20,h=10)
            
    }
    #}}
        
    pDT.ntc_sum_ex_med_diff<-pDT.ntc_sum_ex_med_diff[match(diag_correlation_selected$gene,pDT.ntc_sum_ex_med_diff$gene)]
    
        
    pDT.ntc_sum_ex_med_diff_plot<-merge(pDT.ntc_sum_ex_med_diff,diag_correlation_selected,by="gene")%>% replace(is.na(.), 0)
    
    pDT.ntc_sum_ex_med_diff_plot[,color:="black"]
    pDT.ntc_sum_ex_med_diff_plot <- within(pDT.ntc_sum_ex_med_diff_plot, color[color == "black" & diff == 0] <- 'red')
    
        
    
    g2<-ggplot(pDT.ntc_sum_ex_med_diff_plot, aes(x=ex.vivo, y=reorder(gene,correlation))) + 
      geom_point(size=1, shape=21, fill=pDT.ntc_sum_ex_med_diff_plot$color)+
      geom_segment(data=pDT.ntc_sum_ex_med_diff_plot, mapping=aes(x=ex.vivo, y=reorder(gene,correlation), 
                                                          xend=ex.vivo+diff, yend=reorder(gene,correlation)),arrow = arrow(length = unit(0.09, "cm")), 
                   size=0.09, color="black")+xlim(-2,2)+labs(x="Trajectories(exvivo-invivo)")+FontSize(y.text = 10)+
      theme_bw(12)+theme(axis.title.y = element_blank())
    
    ggsave(out_ex_in("traj_norm_ntc_ex_to_in_arrow.png"),w=20,h=10)
    
    
    source("Ar03_SCRNA_21_01_02_cluster_enrichment_correlation.R")  
    
    diag_correlation_cluster["gene"]<-rownames(diag_correlation_cluster)
    colnames(diag_correlation_cluster)<-c("corr_cluster","gene")
    
    corr_cluster_df<-merge(diag_correlation_selected,diag_correlation_cluster,by="gene",all.x=T)
    corr_cluster_df <- corr_cluster_df %>% replace(is.na(.), 0)
    
 
    
    g3<-ggplot(corr_cluster_df, aes(x=reorder(gene,correlation), y=corr_cluster)) +
      geom_bar(stat='identity') + ylim(-1.000,1.000)+labs(y="Cluster-Correlation")+
      theme_bw(12)+FontSize(y.text = 10)+
      coord_flip()+theme(axis.title.y = element_blank())
    
    ggsave(out_ex_in("correlation_cluster.png"))
    
    g1+g2+g3+ plot_annotation(
      title = 'Exvivo-invivo')
    ggsave(out_ex_in("correlation_cluster_and_traj_norm_ntc_ex_to_in_arrow1.png"))
    
    # logfc correlation and trajectories
    
    pDT.ntc_sum_ex<-pDT.ntc_sum_ex[,c("gene","tissue","median_ntc_norm_traj")]
    
    pDT.ntc_sum_ex_plot<-dcast.data.table(pDT.ntc_sum_ex, gene ~ tissue, value.var = "median_ntc_norm_traj")
    
    pDT.ntc_sum_ex_plot_color<-merge(pDT.ntc_sum_ex_plot,diag_correlation_selected,by="gene")
    
    NTC<-pDT.ntc_sum_ex_plot[gene=="NTC"]
    correlation<-0.00
    NTC<-cbind(NTC,correlation)
    
    pDT.ntc_sum_ex_plot_color<-rbind(pDT.ntc_sum_ex_plot_color,NTC)
    
    #scatter_plot
    
    pDTA<-pDT.ntc_sum_ex_plot_color
    
    
    
    #color--NTC
    ggplot(pDTA, aes(x=ex.vivo, y=in.vivo)) + 
      geom_point(aes(color = pDTA$correlation))+
      
      theme_bw(12) +
      geom_text_repel(aes(label=gene))
    ggsave(out_ex_in("Scatter_Mye_ex_in.vivo_traj_ntc_mean_sd_norm_color_corrlogfc.pdf"),w=8,h=8)
    
    #leuk_ex
    
    pDT.ntc_leuk_sum_med_diff<- pDT.ntc_leuk_sum_med_diff[match(diag_correlation_leuk_selected$gene, pDT.ntc_leuk_sum_med_diff$gene)]
    
    head( pDT.ntc_leuk_sum_med_diff)
    head(diag_correlation_leuk_selected)
    
    pDT.ntc_leuk_sum_med_diff_plot<-merge( pDT.ntc_leuk_sum_med_diff,diag_correlation_leuk_selected,by="gene")%>% replace(is.na(.), 0)
    
    pDT.ntc_leuk_sum_med_diff_plot[,color:="black"]
    
    pDT.ntc_leuk_sum_med_diff_plot <- within(pDT.ntc_leuk_sum_med_diff_plot, color[color == "black" & diff == 0] <- 'red')
    g21<-ggplot( pDT.ntc_leuk_sum_med_diff_plot, aes(x=ex.vivo, y=reorder(gene,correlation))) + 
      geom_point(size=1, shape=21, fill=pDT.ntc_leuk_sum_med_diff_plot$color)+
      geom_segment(data= pDT.ntc_leuk_sum_med_diff_plot, mapping=aes(x=ex.vivo, y=reorder(gene,correlation), 
                                                             xend=ex.vivo+diff, yend=reorder(gene,correlation)),arrow = arrow(length = unit(0.09, "cm")), 
                   size=0.09, color="black")+xlim(-2,2)+labs(x="Trajectories(ex-leuk)")+FontSize(y.text = 10)+
      theme_bw(12)+theme(axis.title.y = element_blank())
      
    ggsave(out_ex_leuk("traj_norm_ntc_ex_to_leuk_arrow.png"),w=8,h=8)

    diag_correlation_leuk_cluster["gene"]<-rownames(diag_correlation_leuk_cluster)
    colnames(diag_correlation_leuk_cluster)<-c("corr_cluster","gene")

    corr_cluster_df<-merge(diag_correlation_leuk_selected,diag_correlation_leuk_cluster,by="gene",all.x=T)
    corr_cluster_df <- corr_cluster_df %>% replace(is.na(.), 0)

   


    g31<-ggplot(corr[corr$gene %in% genes,], aes(x=reorder(gene,correlation), y=corr_cluster)) +
    geom_bar(stat='identity') + ylim(-1.000,1.000)+
    coord_flip()+labs(y="Correlation-clusters")+FontSize(y.text = 10)+
    theme_bw(12)+theme(axis.title.y = element_blank())
    
    t3<-merge(corr,diag_correlation_selected,by="gene")
    colnames(t3)<-c("gene","correlation","corr_cluster","corr_inv")
    head(t3)
    g31<-ggplot(t3, aes(x=reorder(gene,corr_inv), y=corr_cluster)) +
      geom_bar(stat='identity') + ylim(-1.000,1.000)+
      coord_flip()+labs(y="Correlation-clusters")+FontSize(y.text = 10)+
      theme_bw(12)+theme(axis.title.y = element_blank())

    g11+g21+g31+ plot_annotation(
    title = 'Exvivo-leukemia')
    ggsave(out_ex_leuk("correlation_cluster_and_traj_norm_ntc_ex_to_in_arrow1_reordered.png"))

###########################################################
