
add_outcome_variable <- function(read_counts, sample_names) {
  setDT(read_counts)

  read_counts[,
    outcome := do.call(paste, c(.SD, sep = ",")),
    .SDcols = sample_names
  ]

  read_counts[order(outcome)]

  return(read_counts)
}

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


join_merge_data <- function(read_counts) {
  read_counts <-
    read_counts %>%
    map(setDT) %>%
    reduce(
      merge,
      all = TRUE,
      sort = FALSE,
      no.dups = TRUE,
      by = c("cell", "gene", "umi")
    ) %>%
    replace(is.na(.), 0)
  return(read_counts)
}

discard_low_umi_cells <- function(read_counts_dt, min_umi_cell){

  cells_keep <-
    read_counts_dt %>%
    group_by(cell) %>%
    summarize(umis=n()) %>%
    filter(umis >= min_umi_cell) %>%
    pull(cell)

  read_counts_dt <-
    read_counts_dt %>%
    filter(cell %in% cells_keep)
  return(read_counts_dt)
}


create_datatable <- function(out, min_umi_cell) {
  nsamples <- length(out$sample_names)
  read_counts <- list()
  for (i in seq_len(nsamples)) {
    read_counts[[i]] <- tibble(
      cell = out$cells[[i]],
      umi = out$umis[[i]],
      gene = out$genes[[i]],
      reads = out$nreads[[i]]
    )
    if(min_umi_cell> 1) {
    read_counts[[i]] <-
      discard_low_umi_cells(read_counts[[i]],
                            min_umi_cell)
    }
  }
  names(read_counts) <- out$sample_names
  read_counts <- rename_var_data_list(read_counts, out$sample_names)
  return(read_counts)
}

#' Step 3: Create a datatable of gene-mapped read counts
#' @param out A list of read counts data for all samples.
#' @param min_umi_cell The minimum number of UMIs associated with a cell barcode. Barcodes with fewer UMIs are discarded.
#'  Default is 1 (i.e. all cell barcodes are retained).
#'  For large datasets with many samples set to a larger value (50-200) to reduce memory requirements.
#' @return A list *out* containing sample_names and read_counts (a datatable of joined read counts table with an outcome variable column)
#' @details The read counts is created by joining data from all samples into a single datatable
#' in which each row contains the mapped read counts across all the samples for each unique cell
#' barcode-gene-umi combination that make up the three first columns. An outcome character
#' variable (e.g. "(3,0,2,1)") to group the observations into unique outcomes is added as the last column.
#' @order 3
#' @export
join_data <- function(out, min_umi_cell=1) {
  tic("Joining data. Step 1: creating read counts datatables from lists")
  read_counts <- create_datatable(out, min_umi_cell)
  toc()
  tic("Joining data. Step 2: joining and merging datatables for all samples keyed by cell, umi, and gene")
  out[1:4] <- NULL
  read_counts <- join_merge_data(read_counts)
  toc()
  tic("Joining data. Step 3: creating an outcome variable for each row")
  read_counts <- add_outcome_variable(read_counts, out$sample_names)
  out$read_counts <- read_counts
  toc()
  return(out)
}
