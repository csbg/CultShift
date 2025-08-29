source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
require(ggrepel)
require(WriteXLS)
require(patchwork)


# FUNCTIONS ---------------------------------------------------------------
ds <- function(path){load(path); return(monocle.obj)}

# SETTINGS ----------------------------------------------------------------
ff <- list.files(dirout_load("SCRNA_20_Summary")(""))
(TISSUES <- gsub("_.+", "", ff[grepl("_monocle.singleR$", ff)]))
SIGS.USE <- fread("metadata/markers.signatures.use.scRNA.tsv")
SIGS.USE[, sig := paste(DB, Celltype)]
SIGS.USE.DA <- fread("metadata/markers.signatures.use.scRNA2.tsv")
SIGS.USE.DA[, sig := paste(DB, Celltype)]
# Load data sets
#directories
base <- "FIG_cluster_enrichment/"
basedir <- dirout("FIG_cluster_enrichment/")
InDir5 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
InDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
limmaRes <- read_rds(InDir5("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
meta <- read_rds(InDir5("meta.rds"))
meta$sample1 <- rownames(meta)

# Check if there are at least 3 distinct samples per tissue for each genotype and celltype
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)

selected_KOs <- ko_flags %>%
  group_by(celltype)%>%
  filter(valid_ko) %>%
  pull(genotype) %>% unique()
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

koi <- selected_KOs#coefficients,
###########
# Load annotation ---------------------------------------------------------
inDir <- dirout("Ag_ScRNA_17_proj_celltype_in_monocle_obj")
annList <- fread(inDir("exvivo_Annotations.tsv"))
annList <- rbind(fread(inDir("invivo_Annotations.tsv")),annList)

annList <- annList[!is.na(mixscape_class)]%>%
  filter(celltype_projection %in% c(
    "Gran. P", "GMP", "Gran.", "HSC", "MkP", "Mono", "MEP (early)",
    "Eo/Ba")
  )


# Cell types --------------------------------------------------------------
celltype.ann <- CLEAN.CELLTYPES

# numeric clusters
ann.numeric <- readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle_Clusters.RDS"))
annList$Cluster.number <- ann.numeric[match(annList$rn, rn)]$functional.cluster

# . Cluster enrichment analyses ---------------------------------------------
tx <- "in.vivo"
inDir <- dirout_load("Ag_SCRNA_21_02_ClusterEnrichments_simple")

out <- dirout(paste0(base, "/", "cluster.enrichments/"))
  
# Collect enrichment scores
typex <- "detailed_ex.vivo_withMixscape"
for(typex in gsub("Guides_Fisher_Mixscape_(.+).pdf",
                  "\\1", list.files(inDir(""),
                                    pattern="Guides_Fisher_Mixscape_.*.pdf"))){
    fish.file <- inDir("Guides_Fisher_Mixscape_",typex,".tsv")
    if(!file.exists(fish.file)) next
    
    fish.full <- fread(fish.file)
    fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]
    grep("S", unique(fish.full$mixscape_class), value = TRUE)
    #fish.full <- merge(fish.full, unique(SANN[,c("sample_broad", "timepoint"),with=F]), by.x="sample", by.y="sample_broad")
    timex <- "9d"
    for(timex in c(unique(fish.full$sample))){
      fish <- copy(fish.full)
      #fish <- fish[mixscape_class == "Rbbp4"]
      if(timex != "all") fish <- fish[sample == timex]
      
      # summarize across NTCs
      fish <- fish[, .(
        log2OR = mean(log2OR), 
        dir = length(unique(sign(log2OR[padj < 0.01]))) <= 1, 
        #dir=length(unique(sign(log2OR)))==1, 
        padj = sum(padj < 0.01),
        N =.N,
        Ncells = sum(unique(guide.cells))), by=c("sample", "Clusters", "mixscape_class")]
      fish[dir == FALSE, padj := 0]
      fish[dir == FALSE, log2OR := NA]
      
      # legacy
      fish[, gene := mixscape_class]
   
      fish[padj == 0 | is.na(log2OR), log2OR := 0]
      fish[, sig.frac := padj / N]
      fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
      fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "Clusters", value.var = "log2OR")
      fish[, Clusters := gsub("^Gran$", "Gran.", Clusters)]
      #fish[, Clusters := cleanCelltypes(Clusters)]
      #fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")
      ggplot(fish, aes(x=gene, y=Clusters, size=sig.frac, color=log2OR_cap)) + 
        themeNF(rotate=TRUE) +
        scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
        scale_size_continuous(name="% sign.", range = c(0,5)) +
        geom_point() +
        geom_point(shape=1, color="lightgrey") +
        xlab("Gene") + ylab("Cell type")
      ggsaveNF(
        out("Cluster_enrichments_",typex,"_", timex, ".pdf"), 
        w=length(unique(fish$gene))*0.05 + 0.5,
        h=length(unique(fish$Clusters))*0.05 + 0.5)
      write.tsv(fish, out("Cluster_enrichments_",typex,"_", timex,".tsv"))
    }
  }
#}

# Identify `withMixscape` datasets
withMixscape_files <- grep("_withMixscape", 
                           list.files(inDir(""),
   pattern="Guides_Fisher_Mixscape_detailed.*.withMixscape.tsv"),
                           value = T)

# Read all `withMixscape` data and combine

fish.full <- rbindlist(lapply(withMixscape_files, function(file) {
  dt <- fread(inDir(file))  # Read file
  dt[, tissue := gsub("Guides_Fisher_Mixscape_detailed_(.+)_withMixscape.tsv", "\\1", file)]  # Extract tissue
  return(dt)
}), fill=TRUE)


# Standardize Pu.1 to Spi1
fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]

# Summarize across all `withMixscape` samples

fish <- fish.full[, .(
  log2OR = mean(log2OR, na.rm=TRUE), 
  dir = length(unique(sign(log2OR[padj < 0.01]))) <= 1,  
  padj = sum(padj < 0.01),  # Count significant p-values
  N = .N,
  Ncells = sum(unique(guide.cells))
), by = c("Clusters", "mixscape_class", "tissue")]

# Adjust log2OR for visualization
fish[dir == FALSE, padj := 0]
fish[dir == FALSE, log2OR := NA]

# Setup for plotting
fish[padj == 0 | is.na(log2OR), log2OR := 0]
fish[, sig.frac := padj / N]  # Proportion of significant cases
fish[, log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
fish[, order_key := paste(tissue, Clusters,  sep = "_")]
fish <- hierarch.ordering(dt = fish, toOrder = "mixscape_class", orderBy = "order_key", value.var = "log2OR")

#fish <- hierarch.ordering(dt = fish, toOrder = "mixscape_class", orderBy = c("Clusters","tissue","sample"), value.var = "log2OR")
fish <- fish %>%
  filter(mixscape_class %in% koi)
mapping <- c("GMP" = "GMP", 
             "HSC" = "HSC", 
             "Eo/Ba" = "Eo.Ba", 
             "MkP" = "MkP", 
             "MEP (early)" = "MEP.early", 
             "Mono" = "Mono", 
             "Gran." = "Gran.", 
             "Gran. P" = "Gran.P")

# # Rename 'Clusters' based on the mapping
fish <- fish %>%
   mutate(Clusters = recode(Clusters, !!!mapping))
# 
# # Count occurrences of ex vivo and in vivo for each genotype
# genotype_counts <- fish %>%
#   group_by(mixscape_class, tissue) %>%
#   summarise(n = n(), log2OR_positive = sum(log2OR > 0), .groups = "drop") %>%
#   pivot_wider(names_from = tissue, values_from = c(n, log2OR_positive), values_fill = 0)
# 
# # Filter genotypes where:
# # 1. ex vivo AND in vivo occur >1 times
# # 2. At least one ex vivo and one in vivo have log2OR > 0
# valid_genotypes <- genotype_counts %>%
#   filter(n_ex.vivo > 1 & n_in.vivo > 1 & log2OR_positive_ex.vivo > 0 & log2OR_positive_in.vivo > 0) %>%
#   pull(mixscape_class)
# 
# # Keep only the valid genotypes in the original dataset
# fish <- fish %>%
#   filter(mixscape_class %in% valid_genotypes)

# Check the renamed clusters
unique(fish$mixscape_class)

fish %>% write_rds(basedir("detailed_filtered_cluster_enrichment.rds")) 
# Plot summary
ggplot(fish, aes(y=Clusters, x=tissue, size=sig.frac, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(name="log2(OR)", low="blue", midpoint=0, high="red") +
  scale_size_continuous(name="sign.frac.", range=c(0,5)) +
  geom_point() +
  facet_grid(cols = vars(mixscape_class))+
  geom_point(shape=1, color="lightgrey") +
  xlab("Gene") + ylab("Cell type")+
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0))

############################################
# #broad
# # Identify `withMixscape` datasets
# withMixscape_files <- grep("_withMixscape", 
#                            list.files(inDir(""),
#                                       pattern="Guides_Fisher_Mixscape_broad.*.withMixscape.tsv"),
#                            value = T)
# 
# # Read all `withMixscape` data and combine broad
# 
# fish.full <- rbindlist(lapply(withMixscape_files, function(file) {
#   dt <- fread(inDir(file))  # Read file
#   dt[, tissue := gsub("Guides_Fisher_Mixscape_broad_(.+)_withMixscape.tsv", "\\1", file)]  # Extract tissue
#   return(dt)
# }), fill=TRUE)
# 
# 
# # Standardize Pu.1 to Spi1
# fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]
# 
# # Summarize across all `withMixscape` samples
# 
# fish <- fish.full[, .(
#   log2OR = mean(log2OR, na.rm=TRUE), 
#   dir = length(unique(sign(log2OR[padj < 0.01]))) <= 1,  
#   padj = sum(padj < 0.01),  # Count significant p-values
#   N = .N,
#   Ncells = sum(unique(guide.cells))
# ), by = c("Clusters", "mixscape_class", "tissue")]
# 
# # Adjust log2OR for visualization
# fish[dir == FALSE, padj := 0]
# fish[dir == FALSE, log2OR := NA]
# 
# # Setup for plotting
# fish[padj == 0 | is.na(log2OR), log2OR := 0]
# fish[, sig.frac := padj / N]  # Proportion of significant cases
# fish[, log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
# fish[, order_key := paste(tissue,  Clusters,  sep = "_")]
# fish <- hierarch.ordering(dt = fish, toOrder = "mixscape_class", orderBy = "order_key", value.var = "log2OR")
# 
# #fish <- hierarch.ordering(dt = fish, toOrder = "mixscape_class", orderBy = c("Clusters","tissue","sample"), value.var = "log2OR")
# fish <- fish %>%
#   filter(mixscape_class %in% koi)
# # Plot summary
# ggplot(fish, aes(y=Clusters, x=tissue, size=sig.frac, color=log2OR_cap)) + 
#   themeNF(rotate=TRUE) +
#   scale_color_gradient2(name="log2(OR)", low="blue", midpoint=0, high="red") +
#   scale_size_continuous(name="% sign.", range=c(0,5)) +
#   geom_point() +
#   facet_grid(cols = vars(mixscape_class))+
#   geom_point(shape=1, color="lightgrey") +
#   xlab("Gene") + ylab("Cell type")+
#   optimized_theme_fig()+
#   theme(strip.text.y = element_text(angle = 0))
# 
