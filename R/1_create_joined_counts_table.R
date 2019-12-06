

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


#' Rename variables
#' @param read_counts A list of read counts of all the samples
#' @return A joined read counts table
rename_var_data_list <- function(read_counts) {
  sample_names <- names(read_counts)

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


#' Add outcome variable
#' @param read_counts A list of read counts of all the samples
#' @return A dataframe with an outcome variable added
#' @export
add_outcome_variable <- function(read_counts) {

  S <- ncol(read_counts) - 3
  sample_names <- colnames(read_counts)[4:(S + 3)]

  setDT(read_counts)

  read_counts[,
    outcome := do.call(paste, c(.SD, sep = ",")),
    .SDcols = sample_names
  ]

  read_counts[order(outcome)]

  return(read_counts)
}

create_datatable <- function(cells = cells, umis = umis, genes = genes, nreads = nreads){
  nsamples <- length(cells)
  data <- list()
  for (i in seq_len(nsamples)) {
    data[[i]] <- tibble(cell =cells[[i]],
                        umi = umis[[i]],
                        gene=genes[[i]],
                        reads= nreads[[i]])
  }

  names(data) <- names(cells)
  return(data)
}

read10xMolInfoSamples <- function (samples, barcode.length = NULL)
{
  ref.genes <- NULL
  cells <- umis <- genes <- nreads <- vector("list", length(samples))
  names(cells) <- names(samples)
  for (i in seq_along(samples)) {
    mol.info <- DropletUtils::read10xMolInfo(samples[i], barcode.length = barcode.length)
    if (is.null(ref.genes)) {
      ref.genes <- mol.info$genes
    }
    else if (!identical(ref.genes, mol.info$genes)) {
      stop("gene information differs between samples")
    }
    curgems <- mol.info$data$gem_group
    if (length(curgems) && min(curgems) != max(curgems)) {
      warning(paste0("sample '", samples[i], "' contains multiple GEM groups"))
    }
    current <- mol.info$data
    cells[[i]] <- current$cell
    umis[[i]] <- current$umi
    genes[[i]] <- current$gene
    nreads[[i]] <- current$reads
  }
  return(list(cells = cells, umis = umis, genes = genes, nreads = nreads, ref.genes = ref.genes))
}


#' Join transciptome mapped data by cell, umi, and gene
#' @param read_counts A list of read counts of all the samples
#' @param keys keys for joining
#' @return A joined read counts table
join_read_counts <- function(cells = cells, umis = umis, genes = genes, nreads = nreads, keys = c("cell", "gene", "umi")) {

  read_counts <- create_datatable(cells = cells, umis = umis, genes = genes, nreads = nreads)

  read_counts <-
    read_counts %>%
    rename_var_data_list() %>%
    map(setDT) %>%
    reduce(merge,
           all = TRUE,
           sort = FALSE,
           no.dups = TRUE,
           by = keys
    ) %>%
    replace(is.na(.), 0)

  read_counts$gene <- ref.genes[read_counts$gene]

  return(read_counts)
}

#' Join read counts of all samples
#' @param read_counts A list of read counts of all the samples
#' @param joined_counts_filepath If provided, the filepath of the joined readcounts
#' @return A joined dataframe of read counts and an outcome variable
#' @export
get_joined_read_counts <- function(samples, joined_counts_filepath = NULL) {
  if (!is.null(joined_counts_filepath) & file.exists(joined_counts_filepath)) {
    read_counts <- readRDS(joined_counts_filepath)
  } else if ((!is.null(joined_counts_filepath) & !file.exists(joined_counts_filepath)) | is.null(joined_counts_filepath)) {
    c(cells, umis, genes, nreads, ref.genes) %<-% read10xMolInfoSamples(samples)
    read_counts <- join_read_counts(cells, umis, genes, nreads, ref.genes)
    read_counts <- add_outcome_variable(read_counts)
  }

  if (!is.null(joined_counts_filepath) & !file.exists(joined_counts_filepath)) {
    saveRDS(read_counts, joined_counts_filepath)
  }
  return(read_counts)
}
