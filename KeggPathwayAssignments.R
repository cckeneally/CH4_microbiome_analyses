#Load packges
library(httr)
library(readr)
library(parallel)
library(data.table)
library(progress)
library(KEGGREST)

get_kegg_pathways <- function(ortholog_id) {
  if (is.na(ortholog_id) || ortholog_id == "") {
    return(character(0))
  }
  
  response <- try(keggGet(paste0("ko:", ortholog_id)), silent = TRUE)
  
  if (inherits(response, "try-error")) {
    warning("Error fetching pathways for ortholog ", ortholog_id, ": ", response)
    return(character(0))
  }
  
  pathways <- unlist(lapply(response[[1]]$PATHWAY, function(pathway) {
    return(pathway$entry_id)
  }))
  
  # Add print statements for debugging
  cat("Ortholog ID:", ortholog_id, "\n")
  cat("Pathways:", paste(pathways, collapse = ", "), "\n\n")
  
  return(pathways)
}


#Function to assign pathways to table
#Takes output from picrust 2 full pipeline: /KO_metagenome_out/pred_metagenome_unstrat.tsv

assign_kegg_pathways <- function(input_file, output_file) {
  data <- fread(input_file, sep=",")
  
  # Create a progress bar
  pb <- progress_bar$new(format = "[:bar] :percent :elapsed (ETA: :eta)", total = nrow(data), width = 60)
  
  # Define a wrapper function to update the progress bar
  get_kegg_pathways_with_progress <- function(ortholog_id) {
    result <- get_kegg_pathways(ortholog_id)
    pb$tick()
    return(result)
  }
  
  # Use sapply with the wrapper function
  data$KEGG_Pathways <- sapply(data[[1]], get_kegg_pathways_with_progress, USE.NAMES = FALSE) # assuming the ortholog ID is in the first column
  
  fwrite(data, output_file, sep="\t", quote=FALSE)
}

# Replace below with the appropriate file names
setwd("~/Documents/Postgrad/Data/PhoenixOutputs/coorong16sanalysis_2021/exported-feature-table/picrust2_out_pipeline/KO_metagenome_out")
assign_kegg_pathways("pred_metagenome_test.csv", "pred_metagenome_unstrat_assigned1.tsv")
