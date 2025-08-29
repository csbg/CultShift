library(httr)
library(readr)
library(dplyr)
library(igraph)
library(curl)
library(readr)
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
#
InDir_NTC <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")


# 1. Define Mouse STRING URLs
url_links <- "https://stringdb-static.org/download/protein.links.detailed.v12.0/10090.protein.links.detailed.v12.0.txt.gz"
url_alias <- "https://stringdb-static.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz"

# 2. Define destination filenames
dest_links <- "10090.protein.links.detailed.v12.0.txt.gz"
dest_alias <- "10090.protein.aliases.v12.0.txt.gz"

# 3. Download with SSL verification disabled
curl_download(url_links, destfile = dest_links, handle = new_handle(ssl_verifyhost = 0, ssl_verifypeer = 0))
curl_download(url_alias, destfile = dest_alias, handle = new_handle(ssl_verifyhost = 0, ssl_verifypeer = 0))

# 4. Load data
string_raw <- read.table(dest_links, header = TRUE, sep = "", stringsAsFactors = FALSE)

alias_raw <- read_tsv(dest_alias, col_names = c("protein", "alias", "source"))
alias_raw %>%
  filter(source %in% c("Ensembl_gene", "Gene name", "Ensembl_MGI", 
                       "Ensembl_external_synonym_MGI", "UniProt_GN_Name")) %>%
  head(20)

# 5. Process STRING network
string_edges <- string_raw %>%
  filter(as.integer(combined_score) >= 700) %>%
  mutate(
    protein1 = str_remove(protein1, "^10090\\."),
    protein2 = str_remove(protein2, "^10090\\.")
  ) %>%
  dplyr::select(protein1, protein2, combined_score)

# 6. Process alias mapping
gene_map <- alias_raw %>%
  filter(source %in% c("Ensembl_gene", "Ensembl", "Gene name", "UniProt_GN_Name")) %>%
  distinct(protein, alias) %>%
  mutate(protein = str_remove(protein, "^10090\\."))
# 1. Use only rows with gene symbol aliases (i.e., exclude Ensembl-style aliases)
gene_map_clean <- gene_map %>%
  filter(!str_detect(alias, "^ENS"))  # Keep only actual gene symbols like "Hoxb9", "Gnai3", etc.

# 2. Join with STRING edges
string_named <- string_edges %>%
  left_join(gene_map_clean, by = c("protein1" = "protein")) %>%
  dplyr::rename(gene1 = alias) %>%
  left_join(gene_map_clean, by = c("protein2" = "protein")) %>%
  dplyr::rename(gene2 = alias) %>%
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  distinct(gene1, gene2, combined_score)

# 3. Build the graph using gene symbols
g <- graph_from_data_frame(
  string_named %>% select(gene1, gene2, combined_score),
  directed = FALSE
)


E(g)$weight <- string_named$combined_score


#gene selection
ex_in_NTC_per_ct <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
filtered_data <- ex_in_NTC_per_ct %>%
  mutate(Regulation = group) %>%
  filter(group != "n.s")%>%
  filter(abs(logFC) > 2)
affected_genes <- filtered_data$gene
#ko targets
InDir_int <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
ko_genes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))%>%
  filter( coef %in% grep("ex.vivo",coef,value = T))%>%
  mutate(genotype = gsub("ex.vivo","",coef))%>%
  pull(genotype)%>%
  unique()

ko_proteins <- gene_map_clean %>%
  filter(alias %in% ko_genes) %>%
  pull(alias) %>%
  unique()

affected_proteins <- gene_map_clean %>%
  filter(alias %in% affected_genes) %>%
  pull(alias) %>%
  unique()

sum(ko_proteins %in% V(g)$name)          # Should be > 0
sum(affected_proteins %in% V(g)$name)    # Should be > 0


#function
get_network_proximity <- function(ko_gene, affected_genes, graph) {
  if (!(ko_gene %in% V(graph)$name)) {
    return(NA)
  }
  
  affected_in_graph <- unique(affected_genes[affected_genes %in% V(graph)$name])
  
  if (length(affected_in_graph) == 0) {
    return(NA)
  }
  
  dists <- distances(graph,
                     v = ko_gene,
                     to = affected_in_graph,
                     weights = 1 / E(graph)$combined_score)  # adjust if using another weight
  return(mean(dists, na.rm = TRUE))
}
#
network_proximity_scores <- sapply(ko_genes, function(gene) {
  get_network_proximity(gene, affected_genes, g)
})

components <- components(g)
table(components$membership[V(g)$name %in% ko_proteins])
table(components$membership[V(g)$name %in% affected_proteins])

affected_in_ko_comp <- affected_proteins[components$membership[match(affected_proteins, V(g)$name)] == 1]
library(igraph)

get_network_proximity_weighted <- function(ko_gene, affected_genes, graph, epsilon = 1e-6) {
  if (!(ko_gene %in% V(graph)$name)) {
    return(NA)
  }
  
  # Get connected components
  components <- components(graph)
  ko_comp <- components$membership[which(V(graph)$name == ko_gene)]
  
  # Affected genes in same component
  affected_in_comp <- affected_genes[affected_genes %in% V(graph)$name]
  affected_in_comp <- affected_in_comp[components$membership[match(affected_in_comp, V(graph)$name)] == ko_comp]
  
  if (length(affected_in_comp) == 0) {
    return(0)
  }
  
  # Calculate distances to affected genes in same component
  dists <- distances(graph,
                     v = ko_gene,
                     to = affected_in_comp,
                     weights = 1 / E(graph)$combined_score)
  
  mean_dist <- mean(dists, na.rm = TRUE)
  n_connected <- length(affected_in_comp)
  
  score <- n_connected / (mean_dist + epsilon)
  return(score)
}

network_proximity_scores_weighted <- sapply(ko_proteins, function(gene) {
  get_network_proximity_weighted(gene, affected_proteins, g)
})
