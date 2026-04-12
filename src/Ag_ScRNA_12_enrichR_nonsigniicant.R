###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(latex2exp)

#####################################################################
#InDir2 <- dirout("Figure2_Mye")
#InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")

base  <-  "Ag_ScRNA_12_fgsea_overlap/"
basedir  <-  dirout("Ag_ScRNA_12_fgsea_overlap")

########################################################################
limmaRes  <-  read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
# --- 1️⃣ Define non-significant genes per cell type and coefficient ---
# limmaRes must have: ensg, celltype, coef, logFC, adj.P.Val
non_sig_genes <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    non_sig = list(ensg[!(abs(logFC) > 1 & adj.P.Val < 0.05)]),
    .groups = "drop"
  )

# --- 2️⃣ Choose enrichR databases ---
enrichr_dbs <- c(
  "KEGG_2019_Mouse",
  "MSigDB_Hallmark_2020",
  "WikiPathways_2019_Mouse",
  "GO_Biological_Process_2021"
)

# --- 3️⃣ Function to run enrichR ---
run_enrichr <- function(gene_list, dbs = enrichr_dbs) {
  if(length(gene_list) == 0) return(NULL)
  enrichr(gene_list, dbs)
}

# --- 4️⃣ Run enrichment for each non-significant gene set ---
non_sig_enrichment <- non_sig_genes %>%
  mutate(enrichr_res = map(non_sig, ~ run_enrichr(.x, enrichr_dbs)))

# --- 5️⃣ Flatten results to one table per database ---
non_sig_enrichment_flat <- non_sig_enrichment %>%
  mutate(enrichr_flat = map(enrichr_res, ~ {
    if(is.null(.x)) return(NULL)
    bind_rows(
      lapply(names(.x), function(db) {
        df <- .x[[db]]
        df$db <- db
        df
      }),
      .id = "source"
    )
  })) %>%
  dplyr::select(celltype, coef, enrichr_flat) %>%
  unnest(enrichr_flat)

# --- 6️⃣ Save results ---
write_rds(non_sig_enrichment_flat, file = basedir("enrichment_non_sig_genes.rds"))
write.table(non_sig_enrichment_flat, file = basedir("enrichment_non_sig_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
non_sig_enrichment_flat_filtered_terms <- non_sig_enrichment_flat %>%
  filter(Odds.Ratio >10) %>%
  filter(Adjusted.P.value < 0.01) %>%
  pull(Term) %>%
  unique()
grep("interaction",non_sig_enrichment_flat$coef, value = T)
gsea.res <- non_sig_enrichment_flat %>% 
  filter(coef %in% unique(grep("interaction",non_sig_enrichment_flat$coef, value = T)))%>%
    filter(Term %in% non_sig_enrichment_flat_filtered_terms)
dbx <- "MSigDB_Hallmark_2020"
pDT <- gsea.res %>%
  filter(db == "MSigDB_Hallmark_2020")%>%
  arrange(Adjusted.P.value)%>%
  top_n(n = 50)


#pDT$celltype <- factor(pDT$celltype,
# levels =c("HSC","MEP.early","MkP" ,
#           "GMP", "Gran.P", "Gran.", "Mono","Eo.Ba" ))}
# 
# Step 3: Plot with the new pathway order (highest NES first)
Fig1C <- ggplot(pDT, aes(x=coef, y=Term, color = log2(Odds.Ratio), size=pmin(5, -log10(Adjusted.P.value)))) +
  
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log2(Odds.Ratio)"))+
  #name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0,1.8),
    #limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  
  
  #xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(x = NULL,
       y = "Pathways",
       title = "Enriched pathways")+
  
  #coord_flip()+
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
  ),
  legend.position = "right", legend.direction = "vertical",
  legend.justification = "bottom")+
  facet_grid(~celltype)

Fig1C
unique(pDT$coef)
