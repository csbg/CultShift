source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)#
library(dplyr)
library(ggplot2)
library(purrr)


#source("~/code/resources/RFunctions/Basics.R")
base<-"IRF2"
out <- dirout("IRF2")
treatment_dir <- dirout(paste0("IRF2/","Treatment_conditions"))
countdata<-read.delim("/media/AGFORTELNY/PROJECTS/TfCf_AG/IRF2/raw_IRF2_gene.csv",sep = ",",header = T)
rownames(countdata) <- countdata$X
countdata$X<-NULL
nrow(countdata)


# Filtering
d0 <- DGEList(countdata)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
nrow(d)
######
#metadata
# Extract sample names
samples <- colnames(countdata)
treatments <- c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
# Define metadata with correct assignment of treatment_time and IRF genotype
metadata <- data.frame(
  sample = samples,
  treatment_time = factor(
    ifelse(
      grepl("IFNb_4h", samples), "IFNb_4h",
      ifelse(
        grepl("IFNb_24h", samples), "IFNb_24h",
        ifelse(
          grepl("IFNg_4h", samples), "IFNg_4h",
          ifelse(
            grepl("IFNg_24h", samples), "IFNg_24h",
            "Ut"  # Default for untreated samples
          )
        )
      )
    ),
    levels = c("Ut", "IFNb_4h", "IFNb_24h", "IFNg_4h", "IFNg_24h")
  ),
  IRF1 = factor(
    ifelse(grepl("IRF1", samples) | grepl("IRF12", samples), "KO", "WT"),
    levels = c("WT", "KO")
  ),
  IRF2 = factor(
    ifelse(grepl("IRF2", samples) | grepl("IRF12", samples), "KO", "WT"),
    levels = c("WT", "KO")
  )
)

# Ensure Rosa is set as wildtype
metadata$IRF1[grepl("Rosa", metadata$sample)] <- "WT"
metadata$IRF2[grepl("Rosa", metadata$sample)] <- "WT"
rownames(metadata) <- metadata$sample

# Ensure that the columns in countdata match the row names in metadata
countdata <- countdata[, rownames(metadata)]
stopifnot(all(colnames(countdata) == rownames(metadata)))
#
# Create the design matrix
metadata %>% write_rds(treatment_dir("metadata.rds"))
design <- model.matrix(~treatment_time * IRF1 * IRF2, data = metadata)


png(out("voom.png"))
dataVoom<-voom(d, design, plot = TRUE)
dev.off()
dataVoom$E%>% write_rds(treatment_dir("dataVoom.rds"))
# Fit model
colnames(design) <- make.names(colnames(design))
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)
colnames(limmaFit)
limmaRes <- list() # start an empty list
for(coefx in colnames(coef(limmaFit))){ # run a loop for each coef
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coef in a list.
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) %>%
    rownames_to_column("ensg")
}
limmaRes <- bind_rows(limmaRes, .id = "coef")
limmaRes$coef<-gsub("treatment_time","",limmaRes$coef)
limmaRes$coef<-gsub("genotype","",limmaRes$coef)
limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    coef %in% c("IRF1KO", "IRF2KO", "IRF1KO.IRF2KO") ~ paste0(coef,"_Ut"),
    TRUE ~ coef
  ))
unique(colnames(limmaFit$coefficients))
# Define contrasts
# Create contrasts for all time points
contrast.mt <- makeContrasts(
  doubleKO_Ut = IRF1KO + IRF2KO + IRF1KO.IRF2KO,
  IRF1KO_IFNb_4h = IRF1KO + treatment_timeIFNb_4h.IRF1KO,
  IRF2KO_IFNb_4h = IRF2KO + treatment_timeIFNb_4h.IRF2KO,
  doubleKO_IFNb_4h = IRF1KO + IRF2KO + IRF1KO.IRF2KO + 
    treatment_timeIFNb_4h.IRF1KO + 
    treatment_timeIFNb_4h.IRF2KO + 
    treatment_timeIFNb_4h.IRF1KO.IRF2KO,
  IRF1KO_IFNg_4h = IRF1KO + treatment_timeIFNg_4h.IRF1KO,
  IRF2KO_IFNg_4h = IRF2KO + treatment_timeIFNg_4h.IRF2KO,
  doubleKO_IFNg_4h = IRF1KO + IRF2KO + IRF1KO.IRF2KO + 
    treatment_timeIFNg_4h.IRF1KO + 
    treatment_timeIFNg_4h.IRF2KO + 
    treatment_timeIFNg_4h.IRF1KO.IRF2KO,
  IRF1KO_IFNg_24h = IRF1KO + treatment_timeIFNg_24h.IRF1KO,
  IRF2KO_IFNg_24h = IRF2KO + treatment_timeIFNg_24h.IRF2KO,
  doubleKO_IFNg_24h = IRF1KO + IRF2KO + IRF1KO.IRF2KO + 
    treatment_timeIFNg_24h.IRF1KO + 
    treatment_timeIFNg_24h.IRF2KO + 
    treatment_timeIFNg_24h.IRF1KO.IRF2KO,
  IRF1KO_IFNb_24h = IRF1KO + treatment_timeIFNb_24h.IRF1KO,
  IRF2KO_IFNb_24h = IRF2KO + treatment_timeIFNb_24h.IRF2KO,
  doubleKO_IFNb_24h = IRF1KO + IRF2KO + IRF1KO.IRF2KO + 
    treatment_timeIFNb_24h.IRF1KO + 
    treatment_timeIFNb_24h.IRF2KO + 
    treatment_timeIFNb_24h.IRF1KO.IRF2KO,
  IRF2_vs_IRF1_Ut = IRF2KO - IRF1KO,
  IRF2_vs_IRF1_IFNb_4h = IRF2KO + treatment_timeIFNb_4h.IRF2KO - IRF1KO - treatment_timeIFNb_4h.IRF1KO,
  IRF2_vs_IRF1_IFNb_24h = IRF2KO + treatment_timeIFNb_24h.IRF2KO - IRF1KO - treatment_timeIFNb_24h.IRF1KO,
  IRF2_vs_IRF1_IFNg_4h = IRF2KO + treatment_timeIFNg_4h.IRF2KO - IRF1KO - treatment_timeIFNg_4h.IRF1KO,
  IRF2_vs_IRF1_IFNg_24h = IRF2KO + treatment_timeIFNg_24h.IRF2KO - IRF1KO - treatment_timeIFNg_24h.IRF1KO,
  levels = colnames(design)
)

# Fit the contrast
limmaFit.contrast <- contrasts.fit(limmaFit, contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)

# Extract results
limmaRes.contrast <- list()
for(coefx in colnames(contrast.mt)){ # run a loop for each contrast
  limmaRes.contrast[[coefx]] <- topTable(limmaFit.contrast, coef=coefx, number = Inf) %>%
    rownames_to_column("ensg") %>%
    mutate(coef = coefx)
}
limmaRes.contrast <- bind_rows(limmaRes.contrast)

# Combine the contrast results with the original results
limmaRes <- bind_rows(limmaRes, limmaRes.contrast)

# Sanitize column names in the combined results
limmaRes$coef <- gsub("treatment_time", "", limmaRes$coef)
limmaRes$coef <- gsub("genotype", "", limmaRes$coef)
limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    coef %in% c("IRF1KO", "IRF2KO", "IRF1KO.IRF2KO") ~ paste0(coef,"_Ut"),
    TRUE ~ coef
  ))

# View the unique coefficients in the combined results
unique(limmaRes$coef)
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
