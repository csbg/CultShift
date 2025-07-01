source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
basedir <- dirout("FIG_02_scRNA_UMAPs_ar/")

require(ggrepel)
require(WriteXLS)
require(patchwork)
inDir1 <- dirout_load("/SCRNA_10_collect_UMAPs")

xu <- xlab("UMAP 1")
yu <- ylab("UMAP 2")

# SETTINGS ----------------------------------------------------------------

annList <- readRDS(inDir1("ProjVivo_celltypes.RDS"))

# Other projections
umap.proj <- list(
  original=readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjMonocle.RDS")),
  izzo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivo.RDS")),
  in.vivo.X = readRDS(dirout_load("SCRNA_10_collect_UMAPs")("ProjVivoX.RDS"))
)

# SIMPLE SETUP ENDS HERE ---------------------------------------------------------

mobjs <- list()

tissue<-c("ex.vivo","in.vivo")
for(tissuex in tissue){
  (base::load(PATHS$SCRNA$MONOCLE.DIR(paste0(tissuex,"/soupx/"))))
  mobjs[[tissuex]] <- monocle.obj
}

# Remove duplets
for(tissuex in names(mobjs)){
  mobjs[[tissuex]] <- mobjs[[tissuex]][,colnames(mobjs[[tissuex]]) %in% annList$rn]
}

# Initialize the guide and tissue columns in annList
annList$guide <- NA
annList$tissue <- NA

# Match and extract guide and tissue information from `ex.vivo`
match_ex_vivo <- match(annList$rn, colnames(mobjs$ex.vivo))
annList$guide[!is.na(match_ex_vivo)] <- mobjs$ex.vivo@colData$guide[match_ex_vivo[!is.na(match_ex_vivo)]]
annList$tissue[!is.na(match_ex_vivo)] <- "ex.vivo"  # Assign "ex.vivo" for matched rows

# Match and extract guide and tissue information from `in.vivo`
match_in_vivo <- match(annList$rn, colnames(mobjs$in.vivo))
annList$guide[is.na(annList$guide) & !is.na(match_in_vivo)] <- mobjs$in.vivo@colData$guide[match_in_vivo[!is.na(match_in_vivo)]]
annList$tissue[is.na(annList$tissue) & !is.na(match_in_vivo)] <- "in.vivo"  # Assign "in.vivo" for matched rows

# Filter annList to keep only rows where guide is "NTC"
annList <- annList[annList$guide %in% c("NTC","Ash1l") ]
unique(annList$functional.cluster)
annList <- annList %>%
  mutate(
    functional.cluster = case_when(
      functional.cluster %in% c("MEP (G1)", "MEP (pert.)", "MEP (S)") ~ "MEP",  # Grouping all MEP categories to "MEP"
      functional.cluster == "Imm. B-cell" ~ "Imm.B.cell",  # Correct "Imm. B-cell" to "Imm.B.cell"
      is.na(functional.cluster) ~ "Unknown",  # Handle NA values if necessary
      TRUE ~ functional.cluster  # Leave other values unchanged
    )
  )

in.vivo.X <- umap.proj$in.vivo.X
# Filter in.vivo.X based on annList rn to keep only NTC samples

in.vivo.X <- inner_join(in.vivo.X, annList, by = c("rn", "tissue"))

# Generate hexbin object based on filtered in.vivo.X (only NTC samples)
hex.obj <- hexbin::hexbin(x = in.vivo.X$UMAP_1, y = in.vivo.X$UMAP_2, xbins = 100, IDs = TRUE)
in.vivo.X <- cbind(in.vivo.X, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
pDT <- in.vivo.X
pDT <- pDT[, .N, by = c("hex.x", "hex.y", "functional.cluster")]
pDT[, sum := sum(N), by = c("hex.x", "hex.y")]
pDT[, frac := N / sum]

# Filter clusters with significant fractions
pDT <- pDT[frac > 0.25]

# Merge summary data back with the original dataset
merged_data <- inner_join(in.vivo.X, pDT, by = c("hex.x", "hex.y", "functional.cluster"), all.x = TRUE)


# Check the unique values in functional.cluster to ensure it's working

# Generate cluster labels for significant clusters
pDT.labels <- pDT[frac > 0.25, .(hex.x = median(hex.x), hex.y = median(hex.y)), by = c("functional.cluster")]
pDT.labels %>% write_rds(basedir("pDT.labels.rds"))
# color coding
cluster_colors <- c(
  "Mono" = "#E69F00",      # Orange
  "Eo/Ba" = "#56B4E9",     # Sky Blue
  "GMP" = "#009E73",       # Green
  "MEP (early)" = "#F0E442", # Yellow
  "MkP" = "#CC79A7",       # Pink/Purple
  "Gran. P" = "#0072B2",   # Blue
  "Gran." = "#D55E00",     # Reddish Orange
  "HSC" = "#A020F0",       # Purple
  "GMP (early)" = "#999999",  # Light Gray (unchanged)
  "CLP" = "#D9D9D9",       # Light gray for CLP
  "unclear" = "#B0B0B0",    # Gray for unclear
  "Imm. B-cell" = "#8DA0CB", # Soft Blue
  "MEP" = "#D3D3D3",        # Lighter gray for MEP
  "Ery" = "#A9A9A9",        # Slightly darker gray for Ery
  "Imm.B.cell" = "gray"     # Other
)


merged_data %>% write_rds(basedir("Cross_projected_on_in.vivo.rds"))
# Ensure factor ordering for correct label display
merged_data$functional.cluster <- factor(merged_data$functional.cluster, 
                                         levels = names(cluster_colors))

ggplot(merged_data[tissue != "leukemia"], aes(x = UMAP_1, y = UMAP_2)) + 
  
  geom_point(aes(color = functional.cluster), size = 0.00000001 ) + 
  
   geom_text_repel(data = pDT.labels %>%
                    filter(functional.cluster %in% c("Mono", "Eo/Ba", "GMP", "MEP (early)",
                                                     "MkP", "Gran. P", "Gran.", "HSC"
                                                     )),
                  aes(x = hex.x, y = hex.y, label = functional.cluster),
                  size = 1,                  # Adjust text size
                  box.padding = 0.21,         # Distance from points
                  point.padding = 0.21,       # Distance from label anchor
                  segment.color = "black",   # Line color
                  segment.size = 0.004,        # Line thickness
                  force = 10,                # Repelling force
                  max.overlaps = Inf) +
  facet_grid(. ~ tissue) + 
  # Defining color manual scale for clusters
  # Defining color manual scale for clusters
  scale_color_manual(name = "Celltype", values = cluster_colors) + 
  
  # Adjusting the alpha scale for fraction
  #scale_alpha_continuous(name = "Fraction", range = c(0, 1)) +  # To use fraction values for transparency
  
  # Adjusting the theme
  optimized_theme_fig() +
  
  # Positioning legends separately
  # Positioning legends separately
  theme(
    legend.position = "bottom",  # Color legend at the bottom
    legend.box = "horizontal",   # Horizontal alignment of legends
    legend.text = element_text(size = 5),     # Adjust legend text font size
    legend.title = element_blank(), # Remove title for color legend
    legend.spacing = unit(0.5, "cm"),  # Adjust spacing between legends
    legend.key.size = unit(0.5, "lines")  # Adjust size of the legend keys
  ) 
# Correct axis labels (you can define xu and yu separately if needed)
ggsave(outBase("UMAP_InvivoX_NTC.pdf"), w=3.5,h=2.5)


ggplot(merged_data[tissue != "leukemia"], aes(x=UMAP_1, UMAP_2)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..)), bins=100) + 
  scale_fill_gradientn(colors = c("lightgrey", "#1f78b4", "#e31a1c", "#ff7f00")) +
  facet_grid(. ~ tissue) +
  xu + yu+
  optimized_theme_fig()
ggsave(outBase("UMAP_InvivoX_NTC_distribution.pdf"), w=3,h=1.5)
