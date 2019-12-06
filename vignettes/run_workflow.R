#library(rhdf5)
#library(DropletUtils) # install but not load
library(tidyverse)
library(matrixStats)
library(broom)
library(furrr)
library(tictoc)
library(data.table)
library(cowplot)
library(scales)
library(zeallot) # For %<-% that unpacks lists in the Python manner
options(future.fork.enable =TRUE)

plan(multiprocess)
dataset_name <-"hiseq4000"
input_dir <- file.path("/home/rfarouni/Documents/index_hopping/data/hiseq4000/input")
output_dir <- file.path("/home/rfarouni/Documents", "hiseq4000")



joined_counts_filepath <-  file.path(output_dir, sprintf("%s_read_counts.rds",  dataset_name))


samples <- get_h5_filenames(input_dir)

out<- purge_phantoms(samples, torc=3, joined_counts_filepath=joined_counts_filepath)




# #dir.create(output_dir)





tic("Step 7: create umi counts matrices")

toc()



read_counts <-
  read_counts %>%
  select(-c("retain", paste0(sample_names, "_hat")))


data_list <-
  list(
    umi_counts_cell_gene = umi_counts_cell_gene,
    read_counts = read_counts,
    outcome_counts = outcome_counts,
    fit_out = fit_out,
    summary_stats = summary_stats
  )

toc()

#' Load CellRanger 10xv3 data
#' @param sample The name of the filename (i.e. sample)
#' @return A dataframe
#' @importFrom rhdf5 h5read
read10xMolInfo_3 <- function(sample) {
  data <- list()
  all.barcodes <- as.vector(h5read(sample, "/barcodes"))
  all.barcodes <-
    sub("-[0-9]+", "", all.barcodes) # remove GEM group.
  data$cell <-
    all.barcodes[as.vector(h5read(sample, "/barcode_idx")) + 1L]
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/feature_idx")) + 1L
  data$reads <- as.vector(h5read(sample, "/count"))
  gene_name <- h5read(sample, "/features/name")
  gene_ids <- h5read(sample, "/features/id")


  return(data = data)
}

#' Load CellRanger 10xv2 data
#' The function uses cxx_get_cell_barcodes from dropUtils
#' @param sample The name of the filename (i.e. sample)
#' @return A dataframe
#' @importFrom rhdf5 h5read
read10xMolInfo_2 <- function(sample) {
  data <- list()

  data$cell <-
    .Call(
      DropletUtils:::cxx_get_cell_barcodes,
      sample,
      "barcode",
      NULL
    )
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/gene")) + 1L
  data$reads <- as.vector(h5read(sample, "/reads"))
  gene_name <- h5read(sample, "/gene_names")
  gene_ids <- h5read(sample, "/gene_ids")

  nonunique <-
    duplicated(gene_name) |
    duplicated(gene_name, fromLast = TRUE)

  gene_name[nonunique] <-
    str_c(gene_name[nonunique],
          gene_ids[nonunique],
          sep = ":"
    )

  data <- do.call(tibble, data)

  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene = gene_name[gene])

  return(data)
}

#' Load CellRanger data  for all samples
#' @param input_dir The filepath of the directory
#' @return A dataframe of read counts for all samples
#' @importFrom rhdf5 h5ls
#' @export
load_molecule_info_data <- function(input_dir) {
  molecule_info_filepaths <- get_h5_filenames(input_dir)

  available <-
    h5ls(molecule_info_filepaths[[1]],
         recursive = FALSE
    )
  version <- if ("barcode_idx" %in% available$name) {
    "3"
  } else {
    "2"
  }

  if (version == "3") {
    read_counts <- future_map(
      molecule_info_filepaths,
      read10xMolInfo_3
    )
  }
  else {
    read_counts <- map(
      molecule_info_filepaths,
      read10xMolInfo_2
    )
  }

  return(read_counts)
}

#' Purge and save read_counts
#' @param read_counts read counts
#' @param dataset_name dataset name
#' @param sample_names sample names
#' @param output_dir output dir
#' @return
#' @export
purge_and_save_read_counts <- function(read_counts,
                                       dataset_name,
                                       sample_names,
                                       output_dir) {
  read_counts_purged <-
    read_counts %>%
    select(-sample_names) %>%
    filter(retain) %>%
    select(-retain) %>%
    rename_at(c(paste0(sample_names, "_hat")), list(~ str_remove(., "_hat$")))

  purged_data_filepath <-
    file.path(
      output_dir,
      sprintf(
        "%s_read_counts_purged.rds",
        dataset_name
      )
    )

  saveRDS(read_counts_purged, purged_data_filepath)
}

