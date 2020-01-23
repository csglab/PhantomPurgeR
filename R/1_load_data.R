
#' Step 1: return named list of sample filepaths
#' @param input_dir The directory in which the *molecule_info.h5* files are located.
#' @return A named list of samples filepaths.
#' @details  The function read *molecule_info.h5* files produced
#' by 10X Genomics CellRanger software found in user-provided *input_dir*. Either rename the files in the *input_dir* folder
#' as *{sample_name}.h5* before or rename the named list after.
#' @order 1
#' @export
get_h5_filenames <- function(input_dir) {
  input_dir <- normalizePath(input_dir)

  metadata <-
    list.files(
      path = input_dir,
      pattern = "h5",
      recursive = TRUE
    ) %>%
    enframe(name = NULL, value = "sample_name_ext") %>%
    mutate(sample_name = tools::file_path_sans_ext(sample_name_ext)) %>%
    mutate(filepath = file.path(input_dir, sample_name_ext))

  samples_filepaths <- metadata$filepath
  names(samples_filepaths) <- metadata$sample_name

  return(samples_filepaths)
}

#' Step 2: load data from the molecule_info.h5 files
#' @param samples A named list of filepaths.
#' @param barcode_length the length of the cell barcode for v2 data.
#' @return lists of cell, gene, umi, read count, and ref_genes variables.
#' @details For each sample, the function loads from the *molecule_info.h5* file the data corresponding
#' to the following four fields: cell-barcode, gene, umi, and read. The genes field containing
#' the reference gene names is also read from the first sample on the assumption that this field
#' is identical across the samples - as should be the case.
#' @order 2
#' @export
read10xMolInfoSamples <- function(samples, barcode_length = NULL) {
  tic("Loading samples data from molecule_info.h5 files produced by 10X Genomics CellRanger software")

  ref_genes <- NULL
  cells <-
    umis <- genes <- nreads <- vector("list", length(samples))
  sample_names <- names(samples)
  names(cells) <- sample_names
  for (i in seq_along(samples)) {
    message(paste0("Loading sample: ", samples[i]))
    mol_info <-
      read10xMolInfo(samples[i], barcode.length = barcode_length)
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

  toc()

  return(
    list(
      cells = cells,
      umis = umis,
      genes = genes,
      nreads = nreads,
      ref_genes = ref_genes,
      sample_names = sample_names
    )
  )
}
