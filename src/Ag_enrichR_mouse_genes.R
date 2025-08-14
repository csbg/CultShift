#enrichr
ENRICHR <- dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
ENRICHR.DBS <-union(ENRICHR.DBS,
                    c("GO_Biological_Process_2021",
                      "TRRUST_Transcription_Factors_2019",
                      "Reactome_2022",
                      "GO_Molecular_Function_2023",
                      "GO_Biological_Process_2023",
                      "CellMarker_2024"))
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})