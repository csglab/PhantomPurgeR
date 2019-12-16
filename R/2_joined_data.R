#' Add outcome variable
#' @param read_counts A list of read counts of all the samples
#' @return A dataframe with an outcome variable added
add_outcome_variable <- function(read_counts, sample_names) {
  setDT(read_counts)

  read_counts[,
              outcome := do.call(paste, c(.SD, sep = ",")),
              .SDcols = sample_names]

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


create_datatable <- function(out) {
  nsamples <- length(out$sample_names)
  read_counts <- list()
  for (i in seq_len(nsamples)) {
    read_counts[[i]] <- tibble(
      cell = out$cells[[i]],
      umi = out$umis[[i]],
      gene = out$genes[[i]],
      reads = out$nreads[[i]]
    )
  }
  names(read_counts) <- out$sample_names
  read_counts <- rename_var_data_list(read_counts, out$sample_names)
  return(read_counts)
}

#' Join transciptome mapped data by cell, umi, and gene
#' @param A list of read counts data for all samples
#' @return A  joined read counts table with an outcome variable column
#' @export
join_data <- function(out) {
  tic("Joining data. Step 1: creating read counts datatables from lists.")
  read_counts <- create_datatable(out)
  toc()
  tic("Joining data. Step 2: joining and merging datatables for all samples keyed by cell, umi, and gene.")
  out[1:4] <- NULL
  read_counts <- join_merge_data(read_counts)
  toc()
  tic("Joining data. Step 3: creating an outcomes variable for each row.")
  read_counts <- add_outcome_variable(read_counts, out$sample_names)
  out$read_counts <- read_counts
  toc()
  return(out)
}
