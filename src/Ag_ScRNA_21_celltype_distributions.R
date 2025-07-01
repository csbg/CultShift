source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
InDir1 <- dirout("SCRNA_10_collect_UMAPs")
InDir2 <- dirout("/Ag_ScRNA_17_proj_celltype_in_monocle_obj")
InDir3 <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
out <-  dirout("Ag_ScRNA_21_celltype_distributions")

# Input files

in.vivo <- fread(InDir2("invivo_Annotations.tsv"))%>%
  select("tissue","sample_broad","rn","celltype_projection","genotype","mixscape_class.global","Phase")
ex.vivo <- fread(InDir2("exvivo_Annotations.tsv"))%>%
  select("tissue","sample_broad","rn","celltype_projection","genotype","mixscape_class.global","Phase")
KO_genes <- read.delim(InDir3("metadata_guide.tsv"),row.names=1)%>%
  group_by(genotype, tissue, celltype) %>%   # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique() %>% na.omit()

#Combined 
combined_annotation  <- rbind(ex.vivo,in.vivo) %>%
  filter(celltype_projection %in% c("Eo/Ba", "GMP",  "HSC", "MkP",  "Mono", "Gran." )) %>%
  filter( genotype %in% KO_genes)

# Count cell types per sample
cell_counts <- combined_annotation %>%
  group_by(tissue, sample_broad, genotype, celltype_projection) %>%
  summarise(count = n(), .groups = "drop")

cell_proportions <- cell_counts %>%
  group_by(sample_broad, genotype, tissue) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

ggplot(cell_proportions %>% filter(genotype == "NTC"), aes(x = sample_broad, y = proportion, fill = celltype_projection)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(cols= vars(tissue), space = "free", scales = "free") +
  labs(y = "Cell Type Proportion", x = "Sample", fill = "Cell Type") +
  scale_fill_manual(values = c("#963257","#3772A7","#CCA42D","#93A9A5","#CACACA","#E6D5C9"))+
  optimized_theme_fig()
ggsave(out("celltype_distribution.pdf"))


#Cell cycle

# Filter the combined_annotation to include only genotype "NTC"
filtered_data <- combined_annotation %>%
  filter(genotype == "NTC") %>%
  filter(tissue %in% c("ex.vivo", "in.vivo"))

# Plot the distribution of Phase across different tissue types
ggplot(filtered_data, aes(x = Phase, fill = tissue)) +
  geom_bar(position = "dodge", stat = "count") +  # Create bar plot
  facet_wrap(~ tissue) +  # Separate the plots by tissue
  labs(title = "Phase Distribution for Genotype NTC",
       x = "Phase",
       y = "Count",
       fill = "Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

