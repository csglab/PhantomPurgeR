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


  grouping_vars <- tibble(r = r,
                          k_chimera = k_chimera)

  return(grouping_vars)
}

#' Add grouping variables
#' @param outcome_counts outcome dataset
#' @param sample_names sample names
#' @return Datatable of outcome counts
add_vars_to_outcome_counts <-
  function(outcome_counts, sample_names) {
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
create_outcome_counts <- function(read_counts, sample_names) {
  outcome_counts <-
    read_counts %>%
    group_by(outcome) %>%
    add_tally() %>%
    slice(1) %>%
    select(-c(1:3)) %>%
    select(outcome, n, sample_names) %>%
    ungroup()

  outcome_counts <- add_vars_to_outcome_counts(outcome_counts,
                                               sample_names)

  return(outcome_counts)
}



#' Create chimera counts datatable
#' @param outcome_counts outcome dataset
#' @param S number of samples
#' @return Dataframe of chimera counts
create_chimera_counts <- function(outcome_counts, S) {
  # creating a data table of k_chimera counts
  chimera_counts <-
    outcome_counts %>%
    group_by(r, k_chimera) %>%
    tally(n, name = "nn") %>%
    bind_rows(tibble(
      r = rep(100000, S),
      k_chimera = rep(1:S),
      nn = 0
    )) %>%
    ungroup()

  # Converting table from long to wide format
  chimera_counts <-
    chimera_counts %>%
    complete(r,
             k_chimera,
             fill = list(nn = 0)) %>%
    spread(k_chimera, nn) %>%
    filter(r != 100000) %>%
    ungroup()

  # adding tally and proportions variables
  chimera_counts <-
    chimera_counts %>%
    mutate(
      n_cugs = rowSums(.[2:(S + 1)]),
      n_chimeras = rowSums(.[3:(S + 1)]),
      p_chimeras = n_chimeras / n_cugs,
      pcum_cugs = cumsum(n_cugs) / sum(n_cugs)
    )

  return(chimera_counts)
}

#' Estimate hopping rate
#' @param outcome_counts outcome dataset
#' @param max_r Maximum PCR duplication level to consider
#' @return list of two dataframes (glm estimates and chimera counts)
#' @importFrom broom confint_tidy tidy
fit_glm <- function(chimera_counts, S, max_r, conf_level = 0.99) {
  fit_dt <-
    chimera_counts %>%
    filter(r <= max_r & r > 1) %>%
    nest(data = everything()) %>%
    mutate(
      fit = map(data, ~ glm(
        cbind(n_cugs - n_chimeras, n_chimeras) ~ -1 + r,
        data = .x,
        family = binomial(link = log)
      )),
      tidied = map(fit, tidy),
      confint_tidied = map(fit, confint_tidy, conf.level = conf_level),
      #  glanced = map(fit, glance),
      max_r = map(data, ~ as.integer(max(.x$r)))
    )

  glm_estimates <-
    fit_dt %>%
    unnest(c(tidied,
             confint_tidied,
             max_r)) %>%
    select(-std.error,-statistic,-p.value) %>%
    mutate_if(is.double,
              .funs = list(p = ~ exp(.))) %>%
    select(max_r,
           estimate_p,
           conf.low_p,
           conf.high_p) %>%
    rename(phat = estimate_p,
           phat_low = conf.low_p,
           phat_high = conf.high_p) %>%
    mutate(SIHR = 1 - phat,
           SBIHR = (1 - phat) * ((4 * S - 1) / (4 * S - 4)))

  return(glm_estimates)
}

#' Update chimera counts data
#' @param chimera_counts A dataframe of chimera counts
#' @param glm_estimates A dataframe of GLM estimates
#' @return A dataframe of chimera counts
update_chimera_counts <- function(chimera_counts, glm_estimates) {
  chimera_counts <-
    chimera_counts %>%
    mutate(
      phat_chimeras = 1 - glm_estimates$phat ^ r,
      phat_chimeras_low = 1 - glm_estimates$phat_low ^ r,
      phat_chimeras_high = 1 - glm_estimates$phat_high ^ r
    )
  return(chimera_counts)
}

#' Estimate hopping rate
#' @param out list
#' @param max_r Maximum PCR duplication level to consider
#' @return out list containing three additional dataframes (outcome_counts, glm estimates and chimera counts)
#' @export
estimate_hopping_rate <- function(out, max_r = NULL) {
  outcome_counts <-
    create_outcome_counts(out$read_counts, out$sample_names)


  S <- length(out$sample_names)

  chimera_counts <- create_chimera_counts(outcome_counts, S)

  if (is.null(max_r)) {
    max_r <-
      chimera_counts %>%
      pull(r) %>%
      max()
  }

  glm_estimates <-
    fit_glm(chimera_counts,
            S = S,
            max_r)

  chimera_counts <-
    update_chimera_counts(chimera_counts,
                          glm_estimates)
  out <- c(
    out,
    list(
      outcome_counts = outcome_counts,
      glm_estimates = glm_estimates,
      chimera_counts = chimera_counts
    )
  )
  return(out)
}
