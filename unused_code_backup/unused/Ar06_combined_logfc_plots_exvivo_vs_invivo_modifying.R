library(data.table)
library(dplyr)

# Function to process data
process_data <- function(input_prefix, condition, comparison_prefix, comparison_condition, output_filename) {
  # Read input data
  input_data <- fread(dirout(paste0(inp, input_prefix, "_", condition, "_useClusters"))
                      ("DEG_Results_all.tsv")) %>% 
    dplyr::select(c("guide", "gene_id", "estimate"))
  
  # If there is a comparison_condition, filter the data
  if (!is.null(comparison_condition)) {
    comparison_data <- fread(dirout(paste0(inp, comparison_prefix, "_", comparison_condition, "_useClusters"))
                             ("DEG_Results_all.tsv")) %>% 
      dplyr::select(c("guide", "gene_id", "estimate")) %>% 
      filter(guide %in% input_data$guide)
  } else {
    comparison_data <- NULL
  }
  
  # Combine input and comparison data
  combined <- merge(input_data, comparison_data, by = c("guide", "gene_id"))
  colnames(combined) <- c("guide", "geneid",input_prefix
                          , comparison_prefix)
  
  # Create a grouping column
  combined["group"] <- "not.sign"
  
  # Assign groups based on conditions
  combined[which(combined[input_prefix] >= 1 & combined[comparison_prefix] >= 1.0), "group"] <- "i:up"
  combined[which(combined[input_prefix] <= -1 & combined[comparison_prefix] <= -1.0), "group"] <- "h:down"
  combined[which(combined[input_prefix] <= -1 & combined[comparison_prefix] >= 1.0), "group"] <- "g:inv_up.ex_down"
  combined[which(combined[comparison_prefix] <= -1 & combined[input_prefix] >= 1.0), "group"] <- "f:inv_down.ex_up"
  combined[which(combined[comparison_prefix] >= -1 & combined[comparison_prefix] <= 1 &
                   combined[input_prefix] >= 1.0), "group"] <- "e:inv_low.ex_up"
  combined[which(combined[comparison_prefix] >= -1 & combined[comparison_prefix] <= 1 &
                   combined[input_prefix] <= -1.0), "group"] <- "d:inv_low.ex_down"
  combined[which(combined[input_prefix] >= -1 & combined[input_prefix] <= 1 &
                   combined[comparison_prefix] >= 1.0), "group"] <- "c:ex_low.in_up"
  combined[which(combined[input_prefix] >= -1 & combined[input_prefix] <= 1 &
                   combined[comparison_prefix] <= -1.0), "group"] <- "b:ex_low.in_down"
  
  # Write the result to a file
  write.table(combined, out(output_filename))
}

# Define input parameters
input_conditions <- c("ex.vivo", "in.vivo", "leukemia")
output_filenames <- c("combined_input_prefix_logfc.tsv", "combined_comparison_prefix_logfc.tsv", "combined_leukemia_logfc.tsv")
group_names <- c("a:n.s", "a:n.s", "a:n.s")
comparison_conditions <- c(inv, ex, NULL)  # Specify NULL for cases without comparison

# Process each set of data
for (i in seq_along(input_conditions)) {
  process_data(input_conditions[i], conditions[i], comparison_conditions[i], output_filenames[i], group_names[i])
}
