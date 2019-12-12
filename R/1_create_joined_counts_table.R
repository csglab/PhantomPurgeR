

#' Return filepaths of CellRanger molecule_info.h5 files in directory
#' @param input_dir The directory in which the molecule_info.h5 files are located
#' @return A named list of filepaths
#' @export
get_h5_filenames <- function(input_dir) {
  metadata <-
    list.files(
      path = input_dir,
      pattern = "h5",
      recursive = TRUE
    ) %>%
    enframe(name = NULL, value = "sample_name_ext") %>%
    mutate(sample_name = tools::file_path_sans_ext(sample_name_ext)) %>%
    mutate(filepath = file.path(input_dir, sample_name_ext))

  molecule_info_filepaths <- metadata$filepath
  names(molecule_info_filepaths) <- metadata$sample_name

  return(molecule_info_filepaths)
}

#' Load data from molecule_info.h5 files produced by 10X Genomics CellRanger software.
#' @param samples A named list of filepaths
#' @param barcode_length the length of the cell barcode for v2 data
#' @return lists of cell, gene, umi, read count, and ref_genes variables
#' @export
read10xMolInfoSamples <- function (samples, barcode_length = NULL)
{
  ref_genes <- NULL
  cells <- umis <- genes <- nreads <- vector("list", length(samples))
  sample_names <- names(samples)
  names(cells) <- sample_names
  for (i in seq_along(samples)) {
    mol_info <- DropletUtils::read10xMolInfo(samples[i], barcode.length = barcode_length)
    if (is.null(ref_genes)) {
      ref_genes <- mol_info$genes
    }
    else if (!identical(ref_genes, mol_info$genes)) {
      stop("gene information differs between samples")
    }
    curgems <- mol_info$data$gem_group
    if (length(curgems) && min(curgems) != max(curgems)) {
      warning(paste0("sample '", samples[i], "' contains multiple GEM groups"))
    }
    current <- mol_info$data
    cells[[i]] <- current$cell
    umis[[i]] <- current$umi
    genes[[i]] <- current$gene
    nreads[[i]] <- current$reads
  }
  return(list(cells = cells, umis = umis, genes = genes, nreads = nreads, ref_genes = ref_genes, sample_names=sample_names))
}




#' Add outcome variable
#' @param read_counts A list of read counts of all the samples
#' @return A dataframe with an outcome variable added
add_outcome_variable <- function(read_counts, sample_names) {

  setDT(read_counts)

  read_counts[,
    outcome := do.call(paste, c(.SD, sep = ",")),
    .SDcols = sample_names
  ]

  read_counts[order(outcome)]

  return(read_counts)
}

#' Rename variables
#' @param read_counts A list of read counts of all the samples
#' @return A joined read counts table
rename_var_data_list <- function(read_counts, sample_names) {

  name_list <- function(x, y) {
    names(x[[1]]) <- y
    return(x[[1]])
  }

  reads_vars <-
    rep("reads", length(sample_names)) %>%
    set_names(sample_names) %>%
    map(as.list) %>%
    map2(., sample_names, name_list)

  rename_var <- function(x, y) {
    x <-
      x %>%
      rename(!!!y)
    return(x)
  }

  read_counts <-
    read_counts %>%
    map2(., reads_vars, rename_var)
}




join_merge_data<- function(read_counts) {
  read_counts <-
    read_counts %>%
    map(setDT) %>%
    reduce(merge,
           all = TRUE,
           sort = FALSE,
           no.dups = TRUE,
           by = c("cell", "gene", "umi")
    ) %>%
    replace(is.na(.), 0)
  return(read_counts)
}


create_datatable <- function(output){
  nsamples <- length(output$sample_names)
  read_counts <- list()
  for (i in seq_len(nsamples)) {
    read_counts[[i]] <- tibble(cell =output$cells[[i]],
                               umi = output$umis[[i]],
                               gene=output$genes[[i]],
                               reads= output$nreads[[i]])
  }
  names(read_counts) <- output$sample_names
  return(read_counts)
}

#' Join transciptome mapped data by cell, umi, and gene
#' @param A list of read counts data for all samples
#' @return A  joined read counts table with an outcome variable column
join_read_counts <- function(output) {
  read_counts <- create_datatable(output)
  output[1:4] <- NULL
  read_counts <- rename_var_data_list(read_counts, output$sample_names)
  read_counts <- join_merge_data(read_counts)
  read_counts <- add_outcome_variable(read_counts, output$sample_names)
  output$read_counts <- read_counts
  return(output)
}


