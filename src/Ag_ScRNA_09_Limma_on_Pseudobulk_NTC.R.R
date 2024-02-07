##################################################
#***Ag_ScRNA_08_Limma_on_Pseudobulk_NTC.R********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#06-12-23
#Aarathy
#############
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library("dplyr")
library(tidyverse)
inDir<-dirout("/ScRNA_08_Pseudobulk_limma")
out<-dirout("ScRNA_09_Pseudobulk_limma_enrichment")
########################
#load data
########################
#metadata
meta<-read.delim(inDir("metadata.tsv"),row.names=1)
genotype<-"gene"
tissue<-"ex.vivo"
celltype<-"Mono."
number<-data.table()
for(i in unique(meta$celltype)){
  for (gene in unique(meta$genotype)){
    for (tis in unique(meta$tissue)){
      number<-rbind(number,data.table(
      genotype=gene,
      tissue=tis,
      celltype=i,
      sample=nrow(meta[genotype==gene & tissue==tis & celltype==i])))
    }
  }
}

number<-number[genotype%in%unique(meta[tissue=="ex.vivo",]$genotype),]
ggplot(number,aes(x=genotype,y=sample))+
       geom_col()+
  facet_grid(rows = vars(number$tissue),cols = vars(number$celltype))+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
write.tsv(number,out("celltypes_numbers.tsv"))
#meta<-as.data.frame(meta)
#counts
counts <- read.delim(inDir("combined_in_ex_counts.tsv"), row.names = 1)
colnames(counts)<-gsub("Eo.Ba","Eo_Ba",colnames(counts))
#meta data
meta<-meta[meta$celltype!="B-cell" & meta$celltype!="CLP" & meta$celltype!="Ery" & meta$celltype!="Gran." & meta$celltype!="unclear",]

#change rownames to match column names of expression data
rownames(meta)<-gsub("Eo/Ba","Eo_Ba",rownames(meta))
meta$cell<-gsub("Eo/Ba","Eo_Ba",meta$cell)
meta$celltype<-gsub("Eo/Ba","Eo_Ba",meta$celltype)
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
rownames(meta)<-gsub(" ",".",rownames(meta))
rownames(meta)<-gsub("\\(",".",rownames(meta))
rownames(meta)<-gsub("\\)",".",rownames(meta))
rownames(meta)<-gsub("-",".",rownames(meta))

#selecting only NTC
NTC_counts<-counts[,grep("NTC",colnames(counts),value = T)]
NTC_meta<-meta[grep("NTC",rownames(meta),value = T),]
counts<-counts[,rownames(NTC_meta)]

stopifnot(all(colnames(counts)==rownames(NTC_meta)))

#filtering
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

group <- interaction(NTC_meta$celltype, NTC_meta$tissue)

plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~0 + group)
rownames(mm)<-rownames(NTC_meta)
pheatmap(mm)

#Normalization
dataVoom <- voom(d, mm, plot = T)
limmaFit <- lmFit(dataVoom, mm)
limmaFit <- eBayes(limmaFit)

################################################################################################

# List of comparisons
comparisons <- list(
  c("groupEBMP.in.vivo", "groupEBMP.ex.vivo", "EBMP"),
  c("groupEo_Ba.in.vivo", "groupEo_Ba.ex.vivo", "Eo_Ba"),
  c("groupGMP.in.vivo", "groupGMP.ex.vivo", "GMP"),
  c("groupHSC.in.vivo", "groupHSC.ex.vivo", "HSC"),
  c("groupMEP.in.vivo", "groupMEP.ex.vivo", "MEP"),
  c("groupMkP.in.vivo", "groupMkP.ex.vivo", "MkP"),
  c("groupMono.in.vivo", "groupMono.ex.vivo", "Mono")
)

# Initialize a list to store top table results
top_table <- list()

# Perform differential expression analysis for each comparison
for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  id <- comp[3]
  
  # Create contrast matrix
  contrast <- paste0(group1, "-", group2)
  
  # Ensure that contrast inputs are accepted
  if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
    stop("Invalid column names for comparison:", group1, " and ", group2)
  }
  
  # Fit model and perform differential expression analysis
  tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
  tmp <- eBayes(tmp)
  top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
}

# Combine the results into a single data frame

top_table_res <- bind_rows(top_table, .id = 'id')%>%
  mutate(genes=gsub("\\...\\d+","",rownames(top_table_res)))

# Create a column 'group' to group the genes as up, down, or ns
top_table_res$group <- ifelse(top_table_res$logFC >= 1 & 
                                top_table_res$adj.P.Val <= 0.05, "up", 
                              ifelse(top_table_res$logFC <= -1 & 
                                       top_table_res$adj.P.Val <= 0.05, "down", "n.s"))


# Plot and save results volcano plots
for (i in unique(top_table_res$id)){
  p<-ggplot(data=top_table_res[top_table_res$id==i,],
            aes(x=logFC,y=-log10(adj.P.Val),col=group))+
    geom_point()+
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3"))+
    ggtitle(i)
  ggsave(out(paste(i,".pdf")))
}
# to have in the same plot
ggplot(data=top_table_res,aes(x=logFC,y=-log10(adj.P.Val),col=group))+
  geom_point()+
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3"))+
  ggtitle("NTC_differential_expression")+
  facet_grid(rows = vars(top_table_res$id))
ggsave(out("NTC_differential_expression.pdf"))

#########################
# Enrichment analysis
#########################

for(coefx in unique(top_table_res$id)){
  
  # Extract genes of interests (GOI) for a given coefficient (see yesterday's example)
  goi <- top_table_res %>% filter(group=="up" & id==coefx) %>%
    pull(genes)
  
  # Add code here to perform enrichment analysis (see yesterday's example)
  enr.res <- enrichr(goi,databases = c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021"))
  
  # The results will be a list, where each entry is one database. We will combine those into one long table
  enr.res <- bind_rows(enr.res, .id="db")
  
  # Store results in the list
  enr.res.list[[coefx]] <- enr.res
}
enr.res.all <- bind_rows(enr.res.list, .id="coefx")
enr.res.all%>% filter(Odds.Ratio > 10 & 
                        Adjusted.P.value<0.01 & 
                        db=="MSigDB_Hallmark_2020")%>%
                        ggplot(aes(x=coefx,y=Term,size=log2(Odds.Ratio),
                        color=-log10(Adjusted.P.value)))+
                        geom_point()+scale_color_gradient2(high="red", low="blue") 
ggsave(out("NTC_DEG_enrichment.pdf"))
############################
goi.enr <- enr.res.all %>%
  filter(Adjusted.P.value < 0.01 & Odds.Ratio > 10) %>%
  pull(Genes) %>%
  strsplit(";") %>%
  unlist()%>%
  unique()
enriched<-top_table_res %>%
#mutate(gene = gmap[ensg,]$external_gene_name) %>%
  filter(toupper(genes) %in% goi.enr)

set1<-enriched %>% 
  filter(adj.P.Val<0.01)
  
ggplot(set1[200:450,],aes(x=id,y=genes,size=-log10(adj.P.Val),
                    color=logFC))+
  geom_point()+
  scale_color_gradient2(high="red", low="blue")+
  theme(text = element_text(size=10))
###################################################
goi.all<-top_table_res %>%
  filter(adj.P.Val<0.01) %>%
  group_by(id) %>%
  top_n(10,logFC) %>%
  pull(genes)

(p.coef <- top_table_res %>%
   filter(genes %in% goi.all) %>%
   ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
   geom_point() +
   scale_color_gradient2(high="red", low="blue") +
   theme_bw())
#############################################
#
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- NTC_meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))
############################################
t<-top_table_res %>% filter(adj.P.Val<0.01)%>%
  top_n(1,abs(logFC))%>%
  pull(genes)

head(t)
test2<-NTC_meta%>%
  mutate(exp=dataVoom$E["Gbp4",])
ggplot(test2,aes(x=tissue,y=exp))+
  geom_point()+facet_grid(.~celltype)
