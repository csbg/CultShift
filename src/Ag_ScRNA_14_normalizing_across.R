###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
#####################################################################
inDir<-dirout("/Ag_ScRNA_08_Pseudobulk")
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/"
basedir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/")

meta<-read.csv(inDir("metadata.csv"),row.names=1)
meta<-meta[grep("NTC",rownames(meta),value = T),]
meta$genotype<-NULL
meta_izzo<-read.csv(inDir("izzo_metadata.csv"),row.names=1)

meta<-rbind(meta,meta_izzo)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

#only select the genotypes present in both tissue conditions


# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
#############
#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T.cell","IMP2","IMP1","Neu")
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
unique(meta[meta$tissue=="ex.vivo",]$celltype)
unique(meta[meta$tissue=="izzo",]$celltype)

meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta<-meta[!grepl("T.cell",rownames(meta)),]

meta[meta$celltype=="Eo",]$celltype<-gsub("Eo","Eo.Ba",meta[meta$celltype=="Eo",]$celltype)
meta[meta$celltype=="Ba",]$celltype<-gsub("Ba","Eo.Ba",meta[meta$celltype=="Ba",]$celltype)

#counts
counts <- read.csv(inDir("ex_counts.csv"), row.names = 1)
counts_izzo<-read.csv(inDir("izzo_counts.csv"), row.names = 1)
counts_in<-read.csv(inDir("in_counts.csv"), row.names = 1)
counts<-cbind(counts_izzo,counts_in[rownames(counts_izzo),],counts[rownames(counts_izzo),])
counts<-counts[,rownames(meta)]
stopifnot(all(colnames(counts)==rownames(meta)))

counts<-counts%>%na.omit()
################
# edgeR based normalization
d0 <- counts[, rownames(meta)]
d0 <- DGEList(d0)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
model_formula<- "~tissue"
model_formula <- as.formula(model_formula)
# Define the model matrix internally based on the provided formula and metadata
modelMatrix <- model.matrix(model_formula, data = meta)

#A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
#A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
#The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs
# voom
dataVoom <- voom(d, modelMatrix, plot = T)

dat.list <- list()
#interferon
selection<-"interferon_alpha_genes"
goi.all<-interferon_alpha_genes

# Find indices of empty character strings
goi.all <- Filter(nzchar, goi.all)
goi.all<-goi.all[goi.all %in% rownames(dataVoom)]
dat.list <- list()


for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
p.vals <- bind_rows(dat.list, .id="genes")

(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ celltype, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample_incl_izzo.pdf")),
       w = 30, h = 15)
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(rows = vars(celltype), space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 5)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample_cell_wise.pdf")),
       w = 40, h = 15)
unique(p.vals$sample1)
unique(p.vals$sample)
unique(p.vals$cell)
p.vals <- bind_rows(dat.list, .id="genes") %>%
  mutate(cell = as.character(cell))
dat<-ggplot(data=p.vals,aes(x=tissue, y=genes, fill=E)) + 
  geom_tile() +
  facet_grid(. ~ celltype, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dat1<-dat$data
ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
#####################################################################

load(PATHS$RESOURCES$Enrichr.mouse)
dbx<-unique(names(enr.terms))[1]


#pathway genes
inflammatory_response_genes<-enr.terms$MSigDB_Hallmark_2020$`Inflammatory Response`
interferon_alpha_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`
interferon_gamma_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`

