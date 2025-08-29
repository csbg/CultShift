InDir1 <- dirout("Figure2")
InDir2 <- dirout("Figure2_Mye")


InDir_Mye <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide_ex_with_Mye/")
data_mye <- read.delim(InDir_Mye("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
InDir_wo <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
data_wo <- read.delim(InDir_wo("combined_in_ex_counts_guide.tsv"), row.names = 1)

woMye <- read_rds(InDir1("expression.rds"))
Mye <- read_rds(InDir2("expression.rds"))
head(Mye)
head(woMye)
Mye %>% filter(genotype == "Rcor1", tissue == "ex.vivo", celltype == "Eo.Ba")%>%
  pull(samples) %>%
  unique()
woMye %>% filter(genotype == "Rcor1" , tissue == "ex.vivo", celltype == "Eo.Ba")%>%
  pull(samples) %>%
  unique()

Mye %>% filter(genotype == "NTC", tissue == "ex.vivo", celltype == "Eo.Ba")%>%
  pull(samples) %>%
  unique()
woMye %>% filter(genotype == "NTC" , tissue == "ex.vivo", celltype == "Eo.Ba")%>%
  pull(samples) %>%
  unique()
head(data_mye$E)
head(data_wo$E)
data_mye <- data_mye$E
data_wo <- data_wo$E
