library(dplyr)
library(openxlsx)

write_zscore_tables <- function(zscore_data, output_file) {
  
  # Ensure "geneset" column exists
  if (!"geneset" %in% colnames(zscore_data)) {
    stop("The z-score data must include a 'geneset' column.")
  }
  
  # Get list of gene sets present
  gene_sets <- unique(zscore_data$geneset)
  
  # Create a workbook
  wb <- createWorkbook()
  
  # Add one sheet per geneset
  for (gs in gene_sets) {
    sheet_data <- zscore_data %>% filter(geneset == gs)
    
    addWorksheet(wb, gs)
    # Freeze first row
    freezePane(wb, gs, firstRow = TRUE)
    writeData(wb, gs, sheet_data)
  }
  

  
  # Save file
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  message("File written to: ", output_file)
}
