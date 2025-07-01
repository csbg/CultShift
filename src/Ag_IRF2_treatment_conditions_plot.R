# Load necessary libraries
source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)

#########################
#load data
#########################
out <- dirout("IRF2")
treatment_dir <- dirout(paste0("IRF2/","Treatment_conditions"))
limmaRes <- read_rds(out("limmaRes_all.rds"))
dataVoom <-read_rds(treatment_dir("dataVoom.rds"))
metadata<- read_rds(treatment_dir("metadata.rds"))

########################
#functions
########################
optimized_theme <- function() {
  theme_bw() +  # Start with a clean, white background theme
    theme(
      # Text Elements
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Centered, larger plot title
      axis.title = element_text(size = 14, face = "bold"),              # Bold axis titles
      axis.text = element_text(size = 12, color = "black"), 
      axis.text.x = element_text(angle = 45, hjust = 1),
      #legend.position = "right",# Axis text size and color
      legend.title = element_text(size = 13, face = "bold"),            # Bold legend title
      legend.text = element_text(size = 11),                            # Legend text size
      
      # Background and Grids
      panel.grid.major = element_line(color = "grey85", size = 0.5),    # Light grey gridlines
      panel.grid.minor = element_blank(),                               # No minor gridlines
      panel.background = element_rect(fill = "white", color = NA),      # White panel background
      plot.background = element_rect(fill = "white", color = NA),       # White plot background
      
      # Facets
      strip.background = element_blank(),                               # Remove strip background for facets
      strip.text = element_text(size = 13, face = "bold"),              # Bold facet strip text
      
      # Margins and Spacing
      plot.margin = margin(10, 10, 10, 10),                             # Add plot margin for better spacing
      legend.key = element_blank(),                                     # Remove background behind legend items
      legend.position = "right"                                         # Legend positioned on the right
    )
}

####################################################################
#Gene expression and logFC plots------------------------------------
####################################################################
# Load patchwork library for combining plots

unique(limmaRes$coef)
# Function to prepare data for plotting
prepare_data_for_plotting <- function(limmaRes, KO, treatment) {
  interaction_term <- paste0(treatment, ".", KO)
  
  # Filter for the interaction term and arrange by absolute logFC to get the top 20 genes
  top_genes <- limmaRes %>%
    filter(coef == interaction_term) %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)
  
  return(top_genes$ensg)
}

# List of KO conditions and treatments
KOs <- c("IRF1KO", "IRF2KO")
treatments <- c("IFNb_4h", "IFNg_4h", "IFNb_24h", "IFNg_24h")

# Prepare a combined data frame for all plots using purrr::map
# Prepare plot data for KOs and WT
all_plot_data <- map_dfr(KOs, function(KO) {
  map_dfr(treatments, function(treatment) {
    # Get top genes for the current KO and treatment combination
    top_genes <- prepare_data_for_plotting(limmaRes, KO, treatment)
    
    # Generate relevant coefficient names
    interaction_coef <- paste0(treatment, ".", KO)  # Interaction term for KO and treatment
    KO_effect_ut_coef <- paste0(KO, "_Ut")  # Effect of KO in untreated condition
    KO_effect_treatment_coef <- paste0(KO, "_", treatment)  # KO-specific effect with treatment
    WT <- paste0(treatment)  # Wild-type condition (just treatment)
    
    # Combine the coefficients
    relevant_coefs <- c(interaction_coef, KO_effect_ut_coef, KO_effect_treatment_coef, WT)
    
    # Filter limmaRes for relevant coefficients and top genes
    limmaRes %>%
      filter(ensg %in% top_genes & coef %in% relevant_coefs) %>%
      # Assign "WT" when coef is WT, otherwise assign KO
      mutate(KO = if_else(coef == treatment, "WT", KO),
             treatment = treatment)
  })
})

# Plot generation, assuming you're using ggplot2

treat_out <- dirout(paste0("IRF2", "/Treatment_conditions/",ko))

###########
#
for (tr in unique(all_plot_data$treatment)) {
  for (ko in unique(all_plot_data$KO)) {
    
    if (ko == "WT") {
      # WT does not have an interaction term
      interaction_term <- tr
    } else {
      # KO condition has interaction term
      interaction_term <- paste0(tr, ".", ko)
    }
    
    # Filter data for the current treatment and KO (or WT)
    all_data <- all_plot_data %>%
      filter(treatment == tr, KO == ko, coef == interaction_term)
    
    write.table(all_data, treat_out(paste0(tr, ".", ko, "interaction.tsv")), sep = "\t")
    
    # Get top 50 genes by abs(logFC)
    top_50 <- all_data %>%
      filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
      arrange(desc(abs(logFC))) %>%
      slice_head(n = 50) %>%
      pull(ensg)
    
    # Prepare data for plotting logFC
    data_for_plot <- all_plot_data %>%
      filter(ensg %in% top_50, KO == ko, treatment == tr)
    
    # Update factor levels depending on whether it's WT or KO
    if (ko == "WT") {
      data_for_plot$coef <- factor(data_for_plot$coef, levels = c(paste0(tr)))
    } else {
      data_for_plot$coef <- factor(data_for_plot$coef, levels = c(paste0(tr, ".", ko),
                                                                  paste0(ko, "_", tr),
                                                                  paste0(ko, "_", "Ut"),
                                                                  paste0(tr)))
    }
    head(data_for_plot)
    
    # Prepare expression data
    columns <- if (ko == "WT") {
      c(paste0(tr, "_Rosa_", 1:3),
        paste0("Ut", "_Rosa_", 1:3))
    } else {
      c(paste0(tr, "_", ko, "_", 1:3),
        paste0(tr, "_Rosa_", 1:3),
        paste0("Ut", "_", ko, "_", 1:3),
        paste0("Ut", "_Rosa_", 1:3))
    }
    
    columns <- gsub("KO", "", columns)  # Remove KO from the names if needed
    data_exp <- dataVoom[top_50, columns]
    
    dat.list <- list()
    meta <- metadata[columns, ]
    for (gg in top_50) {
      dat.list[[gg]] <- meta %>%
        mutate(E = scale(data_exp[gg, ])) %>%
        rownames_to_column("samples") %>%
        remove_rownames()
    }
    
    # Combine expression data for plotting
    p.vals <- bind_rows(dat.list, .id = "ensg")
    
    # Convert data to a matrix format for clustering
    expression_matrix <- p.vals %>%
      select(ensg, samples, E) %>%
      pivot_wider(names_from = samples, values_from = E) %>%
      column_to_rownames("ensg")
    
    # Compute the distance matrix and perform hierarchical clustering
    dist_matrix <- dist(expression_matrix)  # Euclidean distance by default
    hclust_result <- hclust(dist_matrix)    # Hierarchical clustering
    
    # Reorder ensg based on clustering results
    ordered_ensg <- rownames(expression_matrix)[hclust_result$order]
    
    # Update the ensg factor in the original data frame to reflect this order
    p.vals$ensg <- factor(p.vals$ensg, levels = ordered_ensg)
    
    # Plot expression heatmap using ggplot2
    expression_plot <- ggplot(p.vals, aes(x = samples, y = ensg, fill = E)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red") +
      labs(title = paste0(ko, "-", tr, "-interaction expression heatmap")) +
      optimized_theme()
    
    # Plot logFC
    data_for_plot <- as.data.frame(data_for_plot)
    data_for_plot$ensg <- factor(data_for_plot$ensg, levels = ordered_ensg)
    
    logFC_plot <- ggplot(data_for_plot, aes(x = coef, y = ensg, color = logFC, size = -log10(adj.P.Val))) +
      geom_point() +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size(range = c(2, 8), guide = guide_legend(title = "Adj. P-Value (-log10)")) +
      labs(title = paste("Top 50 Genes for Interaction Effect in Treatment:", tr),
           x = "Coefficient",
           y = "Gene",
           color = "Log Fold Change (logFC)") +
      optimized_theme() +
      scale_x_discrete(labels = if (ko == "WT") {
        c(paste0(tr))
      } else {
        c(paste0("Interaction:", tr, ":", ko),
          paste0(ko, "_", tr),
          paste0(ko, "_", "Ut"),
          paste0(tr))
      })
    
    # Combine the two plots using patchwork
    combined_plot <- logFC_plot + expression_plot
    
    # Define output directory and save the plot
    treat_out <- dirout(paste0("IRF2", "/Treatment_conditions/", ko))
    ggsave(treat_out(paste0(tr, ".", ko, "combined_logFC_expression.pdf")), 
           plot = combined_plot, height = 50 * 0.25 + 1, width = 20)
  }
}

#
##


#############################################
#volcano
#############################################
for (tr in unique(all_plot_data$treatment)) {
  for (ko in unique(all_plot_data$KO)) {
    tr_term <- paste0(ko,"_",tr) 
    ko_term <- paste0(ko,"_","Ut")
    
    # Filter data for the current treatment
    
      # Separate the data for each coefficient
     tr_data <- limmaRes %>%
      filter(coef == tr_term) %>%
      select(ensg,  logFC,  adj.P.Val, group)%>%
       mutate(category = "treatment")
     colnames(tr_data)<- gsub("adj.P.Val", "treatment_adj.P.Val",colnames(tr_data))
     colnames(tr_data)<- gsub("logFC","treatment_logFC",colnames(tr_data))
     colnames(tr_data)<- gsub("group","treatment_group",colnames(tr_data))
     #%>%
     # mutate(paste0("logFC","_",tr) = logFC))
     #        paste0("adj.P.Val","_",tr) = adj.P.Val)))
     ko_data <- limmaRes %>%
      filter(coef == ko_term) %>%
      select(ensg, logFC, adj.P.Val,group)%>%
       mutate(category = "Untreated")
     colnames(ko_data)<- gsub("adj.P.Val", "Ut_adj.P.Val",colnames(ko_data))
     colnames(ko_data)<- gsub("logFC","Ut_logFC",colnames(ko_data))
     colnames(ko_data)<- gsub("group", "Ut_group",colnames(ko_data))
    
    # Merge the data on the gene identifier
    volcano <- inner_join(ko_data,tr_data, by = "ensg")
    
    # Inspect the merged data
    
    
    # Create the group column based on specified conditions
    volcano <- volcano %>%
      mutate(group = case_when(
        # both_up condition
        Ut_group == "up" & treatment_group =="up" ~ "both_up",
        
        Ut_group == "up" & treatment_group =="n.s" &
        Ut_logFC > 2 & treatment_logFC < 1 ~ "Ut_up",
        
        # both_down condition
        Ut_group == "down" & treatment_group =="down" ~ "both_down",
        Ut_group == "down" & treatment_group =="n.s" &
          Ut_logFC < -2 & treatment_logFC < 1 ~ "Ut_down",
        
        # Ut_up_treatment_down condition
        Ut_group == "up" & treatment_group =="down" ~ "Ut_up_treatment_down",
        Ut_group == "n.s" & treatment_group =="up" &
          abs(Ut_logFC) < 1 & treatment_logFC > 2 ~ "treatment_up",
        Ut_group == "n.s" & treatment_group =="down" &
          abs(Ut_logFC) < 1 & treatment_logFC < -2 ~ "treatment_down",
        # Ut_down_treatment_up condition
        Ut_group == "down" & treatment_group =="up" ~ "Ut_down_treatment_up",
        
        # Default case
        TRUE ~ "others"
      ))
    
    volcano <- as.data.frame(volcano)
    
    
        
    top_genes_ut <- volcano %>%
      filter(group != "others") %>%
      filter(group %in% c("Ut_up",
                          "Ut_down"))%>%
      filter(abs(Ut_logFC) > 2, abs(treatment_logFC) < 1)%>%
      arrange(desc(abs(Ut_logFC)))%>%
      slice_head(n=10)
    top_genes_treatment <- volcano %>%
      filter(group != "others") %>% 
      filter(group %in% c("treatment_up","treatment_down")
      )%>% filter(abs(Ut_logFC) < 1,abs(treatment_logFC) > 2)%>%
      arrange(desc(abs(treatment_logFC)))%>%
      slice_head(n=10)
    
    # Inspect the data after adding the group column
    
    #################################################
    # Plot the data
    # Function to calculate label positions and endpoints for lines
    calculate_label_ends <- function(data, x_offset, y_offset) {
      data <- mutate(data,
                     label_x = Ut_logFC + x_offset,
                     label_y = treatment_logFC + y_offset,
                     line_x = ifelse(Ut_logFC < label_x, Ut_logFC + 0.2, Ut_logFC - 0.2),
                     line_y = ifelse(treatment_logFC < label_y, treatment_logFC + 0.2, treatment_logFC - 0.2))
      return(data)
    }
    
    # Calculate label end points for each group
    top_genes_ut <- calculate_label_ends(top_genes_ut, x_offset = 0.5, y_offset = 0.5)
    top_genes_treatment <- calculate_label_ends(top_genes_treatment, x_offset = 0.5, y_offset = -0.5)
    
    # Plot the data
    unique(volcano$group)
    ggplot() +
      # Hexbin plot for the "others" group
      stat_bin_hex(data = filter(volcano, group == "others"), 
                   aes(x = Ut_logFC, y = treatment_logFC, fill = ..count..), 
                   bins = 50, color = NA, alpha = 0.7) +
      scale_fill_gradient(low = "lightgrey",  high= "steelblue", limits = c(1, 5000), name = "Gene Count") +
      
      # Overlay points for the main groups
      geom_point(data = filter(volcano, group != "others"), 
                 aes(x = Ut_logFC, y = treatment_logFC, color = group), 
                 alpha = 0.9, size = 2.5) +
      
      # Manually setting colors for groups
      scale_color_manual(values = c(
        "both_up" = "gray", 
        "both_down" = "gray", 
        "Ut_up" = "#D0154E", 
        "Ut_down" = "#4C889C",
        "treatment_up" = "#D0154E",
        "treatment_down" = "#4C889C",
        "Ut_up_treatment_down" = "blue",
        "Ut_down_treatment_up" = "blue"
       ),
        name = "Group") +
      
      # Add labels for top genes with lines using geom_text_repel
      # Add labels for top genes with lines using geom_text_repel
      geom_text_repel(data = top_genes_ut,
                      aes(label = ensg, x = Ut_logFC, y = treatment_logFC), 
                      size = 3, 
                      box.padding = 0.3, 
                      point.padding = 0.3, 
                      segment.color = 'grey50',
                      nudge_x = -0.9, nudge_y = -0.7,  # Adjust label position
                      force = 5,  # Increase force for stronger repelling
                      min.segment.length = unit(0.2, "cm")) +
      
      geom_text_repel(data = top_genes_treatment,
                      aes(label = ensg, x = Ut_logFC, y = treatment_logFC), 
                      size = 3, 
                      box.padding = 0.3, 
                      point.padding = 0.3, 
                      segment.color = 'grey50',
                      nudge_x = 0.9, nudge_y = -0.7,
                      force = 5,
                      min.segment.length = unit(0.2, "cm")) +
      
      # Add lines from points to labels
      # geom_segment(data = top_genes_Ut,
      #              aes(x = Ut_logFC, y = treatment_logFC, xend = line_x, yend = line_y),
      #              color = "black", alpha = 1, size = 0.5, lineend = "round") +
      # 
      # geom_segment(data = top_genes_Ut,
      #              aes(x = Ut_logFC, y = treatment_logFC, xend = line_x, yend = line_y),
      #              color = "black", alpha = 1, size = 0.5, lineend = "round") +
      
      labs(title = paste0(ko,"_Ut_vs_",ko,"_",tr),
           x = ko_term,
           y = tr_term) +
      
      theme_minimal() +
      
      theme(legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10))
ggsave(treat_out(paste0(tr_term,"_",ko,".pdf")))  
  }
}
############################################
#
#IRF1 vs IRF2
#
############################################
