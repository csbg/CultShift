##################################################
#***Ag_ScRNA_13_invivo_izzo_KO_limma********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#19-03-24
#Aarathy
###########
########################################################################
source("src/00_init.R")
inDir<-dirout("/ScRNA_08_Pseudobulk_limma")
base<-"Ag_ScRNA_13_invivo_izzo_limma_and_enrich.R/"
basedir<-dirout("Ag_ScRNA_13_invivo_izzo_limma_and_enrich.R")

########################################################################
#load data and clean metadata
########################################################################
#metadata
meta<-read.delim(inDir("metadata.tsv"),row.names=1)
meta<-meta%>%filter(meta$tissue=="in.vivo")
#meta1 <- meta1[!(meta1$celltype %in% celltypes_to_exclude), ]
meta<-meta[,c("sample","celltype","tissue")]
NTC_meta<-meta[grep("NTC",rownames(meta),value = T),]
#rownames(NTC_meta)<-gsub("NA","",rownames(NTC_meta))

#izzo
meta1<-read.delim(inDir("izzo_metadata.tsv"),row.names=1)
meta1$sample<-gsub(".*(WT[0-9]+).*", "\\1", rownames(meta1))
meta1$tissue<-"izzo"

meta<-rbind(NTC_meta,meta1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]


# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)\\-]", ".", rownames(meta))
#rownames(meta)<-gsub("..",".",rownames(meta))
# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))
rownames(meta) <- gsub("E/B", "E.B", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta[] <- lapply(meta, gsub, pattern = "E/B", replacement = "E.B.")
meta[] <- lapply(meta, gsub, pattern = "E.B.", replacement = "Eo.Ba")

meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta<-meta[!grepl("T.cell",rownames(meta)),]
meta[grepl("E.B.",rownames(meta)),]
rownames(meta)
##############

#counts

counts_izzo<-read.csv(inDir("izzo_counts.csv"), row.names = 1)
counts <-read.delim(inDir("combined_in_ex_counts.tsv"), row.names = 1)
counts_NTC<-counts[rownames(counts_izzo),rownames(meta[meta$tissue=="in.vivo",])]
counts<-cbind(counts_NTC,counts_izzo)
counts<-counts[,rownames(meta)]%>%na.omit()

stopifnot(all(colnames(counts)==rownames(meta)))
################################################
#factors and levels
################################################
#meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("izzo", "in.vivo"))
################################################
#perform DE independent of celltype
model_formula <- "~tissue"
source("src/Ag_ScRNA_13_invivo_izzo_KO_limma_function.R")
limmaRes <- performDE(meta, counts,model_formula="~tissue",plot_title = "in.vivo-izzo")

#################
