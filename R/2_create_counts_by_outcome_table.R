#' Create grouping variables
#' @param outcome_counts outcome dataset
#' @param sample_names sample names
#' @return dataframe of grouping variables
#' @importFrom matrixStats rowSums2 rowCounts
create_grouping_vars <- function(outcome_counts, sample_names) {
  S <- length(sample_names)

  outcome_counts <-
    outcome_counts %>%
    select(sample_names) %>%
    as.matrix()

  r <- as.integer(rowSums2(outcome_counts))
  k_chimera <- as.integer(S - rowCounts(outcome_counts, value = 0))


  grouping_vars <- tibble(
    r = r,
    k_chimera = k_chimera
  )

  return(grouping_vars)
}

#' Add grouping variables
#' @param outcome_counts outcome dataset
#' @param sample_names sample names
#' @return Datatable of outcome counts
add_vars_to_outcome_counts <- function(outcome_counts, sample_names) {
  grouping_vars <-
    create_grouping_vars(outcome_counts, sample_names)
  outcome_counts <-
    bind_cols(outcome_counts, grouping_vars) %>%
    arrange(r, k_chimera) %>%
    select(outcome, n, r, k_chimera, everything(), sample_names)
  return(outcome_counts)
}

#' Create outcome counts
#' @param read_counts A list of read counts of all the samples
#' @param joined_counts_filepath If provided, the filepath of the joined readcounts to be saved
#' @return Datatable of outcome counts
#' @export
create_outcome_counts <- function(read_counts) {

  sample_names <-
    setdiff(
      colnames(read_counts),
      c("cell", "umi", "gene", "outcome")
    )

  S <- length(sample_names)

  outcome_counts <-
    read_counts %>%
    group_by(outcome) %>%
    add_tally() %>%
    slice(1) %>%
    select(-c(1:3)) %>%
    select(outcome, n, sample_names) %>%
    ungroup()

  outcome_counts <- add_vars_to_outcome_counts(
    outcome_counts,
    sample_names
  )

  return(outcome_counts)
}
