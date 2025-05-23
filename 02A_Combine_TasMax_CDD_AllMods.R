# Created by: Zach Popp
# Date created: 12.07.2023
# Date modified: 03.01.2024
# Version no: v1
# 
# Overview: After running the raster processing for all climate models and 
# scenarios, we merge all of these files into one large dataset. We will use
# the dataset to generate the median for contemporary hazard, projected hazard,
# and change in hazard across models and within SSPs.
#
library(dplyr)

# Set input file and pattern
#
model_outdir <- ""  # Where processed metrics from script 1A are stored (in this application, storage was by model)
joined_outdir <- "" # Where to save result (joined across models/ssp)
pattern_in <- ""    # Setting pattern for which file to read in (to account for version pattern inc v5)

# Generate list of all climate scenarios and ssps to be incorporated
#
files_names <- list.files(path = model_outdir, 
                              pattern = pattern_in, recursive = TRUE, full.names = FALSE)

files_processed <- paste0(model_outdir, "/", files_names)

# Loop through files to read in, clean, and bind across models and scenarios
#
for (file_in in c(files_names)) {
  
  # Print name of file being read in
  #
  cat("Processing ", file_in, "\n")
  
  input <- readRDS(paste0(model_outdir, file_in))
  
  # Remove ID columns
  #
  ID_cols <- names(input)[grepl("ID", names(input))]
  input <- input[!c(names(input) %in% ID_cols)]
  
  # Remove geometry
  #
  input <- as.data.frame(input)
  input$geometry <- NULL

  if (file_in == files_names[1]) {
    model_final <- input
  } else {
    model_final <- rbind(model_final, input)
  }
}

# Clean naming - still an annual measure, but easier to keep shorter name
#
colnames(model_final) <- gsub("_annual$", "", colnames(model_final))

# Save output
#
saveRDS(model_final, paste0(joined_outdir, "all_models.rds"))

