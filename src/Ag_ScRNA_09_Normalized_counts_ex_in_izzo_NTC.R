#limma_NTC
#Comparison:ex.vivo-in.vivo

###############
source("src/00_init.R")
library(edgeR)
library(limma)
##################################################################################
inDir<-dirout("Ag_ScRNA_08_Pseudobulk_limma")
inDir1<-dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
base<-"Ag_ScRNA_09_Normalized_counts_ex_in_izzo_NTC/"
basedir<-dirout("Ag_ScRNA_09_Normalized_counts_ex_in_izzo_NTC")
########################
##################################################################################
#load data
########################
#metadata
meta<-read.delim(inDir1("metadata_in_ex_guide.tsv"),row.names=1)
meta_izzo<-read.csv(inDir1("izzo_metadata.csv"),row.names = 1)
#meta<-meta[meta$tissue=="ex.vivo",]
meta<-meta%>%filter(genotype=="NTC")
meta_izzo$genotype<-"WT"
meta_izzo$guide<-"NTC"
meta_izzo$cell<-NULL

meta$mixscape_global<-NULL
meta<-rbind(meta,meta_izzo)

#meta data
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell","Gran.","E/B","T-cell-Cd3d+")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]


# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-/]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
#rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
#meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")

meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta$tissue <- factor(meta$tissue, levels=c("izzo", "ex.vivo","in.vivo"))
###############
#Save the meta
meta%>%write_rds(basedir("meta_all_NTC.rds"))

##########
counts <- read.delim(inDir("combined_in_ex_counts.tsv"), row.names = 1)
counts<-counts[,grep("NTC",colnames(counts),value = T)]
counts<-counts[,rownames(meta[meta$tissue == "ex.vivo"| meta$tissue=="in.vivo",])]

counts_izzo<-read.csv(inDir1("izzo_counts.csv"),row.names = 1)
# Replace "Eo/Ba" with "Eo.Ba" in row names
colnames(counts_izzo) <- gsub("Eo/Ba", "Eo.Ba", colnames(counts_izzo))
colnames(counts_izzo) <- gsub("[\\ \\(\\)-]", ".", colnames(counts_izzo))
counts_izzo<-counts_izzo[,rownames(meta[meta$tissue=="izzo",])]

result <- merge(counts, counts_izzo, by = "row.names", all = TRUE)%>%na.omit()
counts_all<-result
rownames(counts_all)<-counts_all$Row.names
counts_all$Row.names<-NULL

all(rownames(meta)==colnames(counts_all))
##################################
d0 <- DGEList(counts_all)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group <- interaction(meta$celltype,meta$tissue)
unique(group)

mm <- model.matrix(~0 + group)
rownames(mm) <- rownames(meta)

# Normalization
dataVoom <- voom(d, mm)
dataVoom%>% write_rds(basedir("dataVoom_perCT_NTC_izzo_ex_in.rds"))
# limmaFit <- lmFit(dataVoom, mm)
# limmaFit <- eBayes(limmaFit)
# 
# # Generate comparisons based on tissue types
# tissue_type1<-"ex.vivo"
# tissue_type2<-"in.vivo"
# comparisons <- lapply(unique(NTC_meta$celltype), function(celltype) {
#   comp1 <- paste0("group", celltype, ".",tissue_type1)
#   comp2 <- paste0("group", celltype, ".", tissue_type2)
#   list(comp1, comp2, celltype)
# })
# 
# # Initialize list to store top table results
# top_table <- list()
# 
# # Perform differential expression analysis for each comparison
# for (comp_list in comparisons) {
#   group1 <- comp_list[[1]]
#   group2 <- comp_list[[2]]
#   id <- comp_list[[3]]
#   
#   # Create contrast matrix
#   contrast <- paste0(group1, "-", group2)
#   
#   # Debugging: Print contrast matrix
#   print(contrast)
#   
#   # Ensure that contrast inputs are accepted
#   if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
#     stop("Invalid column names for comparison:", group1, " and ", group2)
#   }
#   
#   # Fit model and perform differential expression analysis
#   tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
#   tmp <- eBayes(tmp)
#   top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
# }
# 
# # Bind results and add NTC_metadata
# ex_in_NTC_per_ct <- bind_rows(top_table, .id = 'celltype') %>%
#   mutate(genes = gsub("\\...\\d+", "", rownames(.)))
# 
# # Create a column 'group' to group the genes as up, down, or ns
# ex_in_NTC_per_ct$group <- ifelse(ex_in_NTC_per_ct$logFC >= 1 & 
#                                 ex_in_NTC_per_ct$adj.P.Val <= 0.05, "up", 
#                               ifelse(ex_in_NTC_per_ct$logFC <= -1 & 
#                                        ex_in_NTC_per_ct$adj.P.Val <= 0.05, "down", "n.s"))
# 
# ex_in_NTC_per_ct %>% write_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
# # Plotting and saving
# ex_in_NTC_per_ct %>%
#   ggplot(aes(x = logFC, y = pmin(10, -log10(adj.P.Val)), col = group)) +
#   geom_point() +
#   scale_color_manual(values = c("#5782A7", "#B9B8B6", "#8A264A")) +
#   ggtitle(paste0(tissue_type1, " ", tissue_type2, " NTC D.E")) +
#   facet_wrap(vars(celltype), scales = "free_x") +
#   theme_bw()
# 
# 
# ggsave(basedir(paste0(tissue_type1, "_", tissue_type2, "_NTC_D.E_percelltype.pdf")), 
#       w=20,h=15)
# ##########################################
