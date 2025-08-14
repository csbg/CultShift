source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")


require(tidyverse)

basedir <- dirout("Fig7_Supplementary/")
InDir <- dirout("Figure2/")
# load, format and merge data ---------------------------------------------------

limmaRes_int <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(genes = ensg)%>%
  mutate(coef = gsub("interaction","",coef))%>%
  dplyr::select(-ensg)
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))


merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y) %>%
  dplyr::select(coef,celltype,genes,
                logFC_KO,
                logFC_NTC,
                adj.P.Val_KO,
                adj.P.Val_NTC)
correlation_results <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
  )
# Step 1: Filter the merged data as done in correlation analysis
data_ntc_interaction <- merged_data %>%
  inner_join(correlation_results, by = c("coef","celltype"))


data_logFC_both_conditions <- readRDS(InDir_cor("in.vivo_ex.vivo_logFC.rds"))


data_logFC_both_conditions <- data_logFC_both_conditions %>% 
  dplyr::rename(gene = ensg, knockout = genotype)

data_ntc_interaction <- data_ntc_interaction %>% 
  dplyr::rename(gene = genes, knockout = coef, logFC_interaction = logFC_KO
  )

data_lFCs <- merge(x = data_logFC_both_conditions, y = data_ntc_interaction,
                   by=c("gene", "knockout", "celltype"))

with(data_logFC_both_conditions, table(celltype, knockout))
with(data_ntc_interaction, table(celltype, knockout))
with(data_lFCs, table(celltype, knockout))

colnames(data_lFCs)

# does the interaction show the difference of ex vivo and in vivo?
ggplot(data_lFCs, aes(x=logFC_ex.vivo -logFC_in.vivo, y =logFC_interaction)) + geom_hex()
ggplot(data_lFCs, aes(x=(logFC_ex.vivo - logFC_in.vivo) - logFC_interaction)) + stat_ecdf()
diff <- table(with(data_lFCs, (logFC_ex.vivo -logFC_in.vivo) - logFC_interaction))


# Plot scatterplot --------------------------------------------------------
ggplot(data_lFCs, aes(x=logFC_ex.vivo, y=logFC_in.vivo)) + 
  geom_hex() +
  facet_wrap(~celltype + knockout)
ggsave(basedir("Scatterplot.pdf"), w=20,h=20)


# Plot correlations -------------------------------------------------------
for(cor_method in c("spearman", "pearson")){
  pDT <- data_lFCs %>%
    group_by(celltype, knockout) %>% 
    summarize(
      cor=cor(logFC_ex.vivo, logFC_in.vivo, method=cor_method), 
      sig_ex = sum(adj.P.Val_in.vivo < 0.05),
      sig_in = sum(adj.P.Val_ex.vivo < 0.05)) %>% 
    mutate(sig_max = pmax(sig_ex, sig_in))
  ggplot(pDT, aes(x=knockout, y=celltype, fill=cor, size=sig_max)) +
    geom_point(shape=21) + 
    scale_fill_gradient2(low="blue", high="red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(basedir(paste0("Correlations_baseline_", cor_method,".pdf")), w=6,h=4)
  write_csv(pDT, basedir(paste0("Correlations_baseline_", cor_method,".csv")))
}



# Predict new values ------------------------------------------------------
dir.create(basedir("predictions"))

models=list(
  exVivoOnly = "0 + logFC_ex.vivo",
  exVivoAndInteraction = "0 + logFC_ex.vivo + logFC_ex.vivo:logFC_NTC",
  interactionNoIntercept = "0 + logFC_ex.vivo*logFC_NTC",
  interactionAndIntercept = "logFC_ex.vivo*logFC_NTC"
)

ctx <- "GMP"
kox <- "Spi1"
mx <- "exVivoOnly"
for(mx in names(models)){
  results_lm <- list()
  for(ctx in c(unique(data_lFCs$celltype), "all")){
    for(kox in c(unique(data_lFCs$knockout), "all")){
      set.seed(342)
      data_lm <- data_lFCs
      if(ctx != "all") data_lm <- filter(data_lm, celltype == ctx)
      if(kox != "all") data_lm <- filter(data_lm, knockout == kox)
      if(nrow(data_lm) == 0) next
      splits <- cut(sample(1:nrow(data_lm)), breaks=5, labels=FALSE)
      for(splitx in unique(splits)){
        data_lm_train <- data_lm[splits != splitx,]
        data_lm_test <- data_lm[splits == splitx,]
        fit <- lm(formula = as.formula(paste("logFC_in.vivo ~", models[[mx]])), data=data_lm_train)
        data_lm_test$logFC_predicted <- predict(fit, data_lm_test)
        data_lm_test <- data_lm_test %>% 
          select(gene, knockout, celltype, logFC_in.vivo, logFC_predicted) %>% 
          mutate(model = mx) %>% 
          mutate(celltype_trained = ctx) %>% 
          mutate(knockout_trained = kox)
        results_lm[[paste(ctx, kox, mx, splitx, sep="_")]] <- data_lm_test
      }
    }
  }
  results_lm <- rbindlist(results_lm)
  results_lm %>% write_csv(basedir(paste("predictions/Predictions", mx, ".csv", sep="_")))
}

ff <- list.files(basedir("predictions/"), full.names = TRUE)

results_lm <- lapply(ff, read_csv)
results_lm <- bind_rows(results_lm)

str(results_lm)
model_names = c(
  "exVivoOnly" = "0 + logFC_ex.vivo",
  "interactionAndIntercept" = "logFC_ex.vivo*logFC_NTC"
)
head(results_lm)
cor_method <- "pearson"
for(cor_method in c("spearman", "pearson")){
  perf <- results_lm %>% 
    mutate(ct_all = ifelse(celltype_trained == "all", "all cts", "individual cts")) %>% 
    mutate(ko_all = ifelse(knockout_trained == "all", "all kos", "individual kos")) %>% 
    group_by(celltype, knockout, ct_all,  ko_all, model) %>% 
    summarize(
      cor=cor(logFC_predicted, logFC_in.vivo, method=cor_method), 
      cor=cor(logFC_predicted, logFC_in.vivo, method=cor_method)
    )
  write_csv(perf, basedir(paste(cor_method,"correlations_pred_act.csv")))
  ggplot(perf, aes(x=knockout, y=celltype, fill=cor)) +
    geom_tile() + 
    scale_fill_gradient2(low="blue", high="red") +
    theme_bw() +
    facet_grid(model ~  ko_all,
               labeller = labeller(model = model_names))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(basedir(paste0("Correlation_predictions_details_", cor_method, ".pdf")), w=12,h=8)
  
  baseline_DT <- read_csv(basedir(paste0("Correlations_baseline_", cor_method,".csv")))
  
  perf <- baseline_DT %>% 
    mutate(ct_all = "0NA") %>% 
    mutate(ko_all = "0NA") %>% 
    mutate(model = "baseline") %>% 
    select(colnames(perf)) %>% 
    rbind(perf)
  
  ggplot(perf, aes(x=model, y=cor)) +
    geom_violin(fill="lightblue", color=NA) +
    geom_boxplot(color="grey", fill=NA, coef=Inf) + 
    geom_jitter(height=0, width=0.2, shape=1) +
    theme_bw() +
    facet_grid(. ~  ko_all, scales = "free", space = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0("Correlation_predictions_violin_", cor_method, ".pdf"), w=8,h=4)
}
colnames(pred_act)


#Sup.Fig.7A

pred_act <- read_csv(basedir("pearson correlations_pred_act.csv"))

deg <- read_rds(InDir_cor(paste0("DEGs_per_tissue.rds")))  

deg <- deg %>%
  filter(condition == "In Vivo")%>%
  dplyr::select("celltype", "genotype","num_degs")%>%
  dplyr::rename(knockout = genotype)
pred_act <- merge(pred_act, deg,by = c("celltype", "knockout"))

combine <- pred_act
unique(combine$model)

model_names = c(
  "exVivoOnly" = "0 + logFC_ex.vivo",
  "interactionAndIntercept" = "logFC_ex.vivo * logFC_NTC"
  #"ex vivo" = "baseline"
)
unique(combine$model)

ko_flag <- ko_flags%>%
  dplyr::rename(knockout = genotype)
combine <- combine %>%
  inner_join(ko_flag, by = c("knockout", "celltype"))
  
combine <- combine %>% filter(valid_ko)
combine <- combine %>%
  filter(model %in% c(
    "exVivoOnly",
    "interactionAndIntercept"  
  ))
column_order <- read_rds(InDir("column_order.rds"))
combine$knockout <- factor(combine$knockout,
                           levels = column_order)
#geom_point

sup.fig.7a <- combine %>%
  filter(valid_ko)%>%
  filter(ct_all != "individual cts") %>%
  ggplot(aes(
    x = knockout,
    y = celltype,
    #size = pmin(3,log10(sig_in)),
    fill = cor  #
  )) +
  geom_point(aes(
    size = pmin(3,log10(num_degs)),
    fill = cor  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = expression("Pearson's correlation")
  ) +
  scale_size_continuous(
    range = c(0,2.5),
    #limits = c(0,2.5),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n)))
  )+
  facet_grid(ko_all ~ model,
             space = "free",
             scales = "free",
             labeller = labeller(model = model_names))+
  labs(x = "KOs",
       y = "Cell type",
       title =  "Correlation of KO-effects (actual versus predicted)") +
  optimized_theme_fig()+
  theme(
    
    legend.position  = "bottom"
  )
ggsave(basedir("Sup.Fig.7a_dot.pdf"),plot = sup.fig.7a, w = 18, h = 10,units = "cm")

sup.fig.7a <- combine %>%
  filter(valid_ko)%>%
  filter(ct_all != "individual cts") %>%
  ggplot(aes(
    x = knockout,
    y = celltype,
    #size = pmin(3,log10(sig_in)),
    fill = cor  #
  )) +
  geom_tile(

  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = expression("Pearson's correlation")
  ) +
  
  facet_grid(ko_all ~ model,
             space = "free",
             scales = "free",
             labeller = labeller(model = model_names))+
  labs(x = "KOs",
       y = "Cell type",
       title =  "Correlation of KO-effects (actual versus predicted)") +
  optimized_theme_fig()+
  theme(
    
    legend.position  = "bottom"
  )
ggsave(basedir("Sup.Fig.7a.pdf"),plot = sup.fig.7a, w = 18, h = 10,units = "cm")

ex_in <- read_rds(InDir("correlation_deg.rds")) 
ex_in <- ex_in %>%
  dplyr::select("genotype", "celltype", "correlation")%>%
  dplyr::rename(knockout = genotype,
                cor =correlation
  )%>%
  mutate(ct_all = "all cts",
         ko_all = "individual kos",
         model = "ex vivo")

ex_in <- merge(ex_in, deg,by = c("celltype", "knockout"))
ex_in <- ex_in %>%
  inner_join(ko_flag,by = c("celltype","knockout"))

ex_in$knockout <- factor(ex_in$knockout,
                           levels = column_order)
ex_in %>%
  filter(valid_ko)%>%
  filter(ct_all != "individual cts") %>%
  ggplot() +
  geom_tile(aes(
    x = knockout,
    y = celltype,
    #size = pmin(3,log10(sig_in)),
    fill = cor  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = expression("Pearson's correlation")
  ) +
  scale_size_continuous(
    range = c(0,2.5),
    #limits = c(0,2.5),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n)))
  )+
  
  labs(x = "KOs",
       y = "Cell type",
       title =  "Correlation of KO-effects (in vivo versus ex vivo)") +
  optimized_theme_fig()+
  theme(
    
    legend.justification = "right"
  )
