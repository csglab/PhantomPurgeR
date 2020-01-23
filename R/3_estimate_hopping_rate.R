
get_prop_fugues <- function(p_cugs, phat) {
  upper_thresh <- 10
  r <- 1:upper_thresh # first ten terms are sufficient
  m <- p_cugs[r]
  y <- rep((1 - phat), upper_thresh)
  g <- sum(exp(r * log(y) + log(m)))

  return(g)
}


estimate_pi_r <- function(nu, phat, S, pi_r_eps = 0.000001) {
  pi_r <- (nu * (S - 1) + (phat - 1)) / (S * phat - 1)
  pi_r <- pmax(pi_r, pi_r_eps)

  return(pi_r)
}

estimate_pi_r_hat_matrix <- function(conditional, sample_names, phat) {
  S <- length(sample_names)

  pi_r_hat <-
    conditional %>%
    mutate_at(vars(sample_names),
      estimate_pi_r,
      phat = phat,
      S = S
    ) %>%
    mutate(sum_p = rowSums(.[sample_names])) %>%
    mutate_at(vars(sample_names), ~ . / sum_p) %>%
    select(-sum_p, -max_hop_rate, -n_cugs, -p_cugs)

  return(pi_r_hat)
}

compute_summary_stats <- function(outcome_counts, phat) {
  sample_names <-
    setdiff(
      colnames(outcome_counts),
      c("n", "r", "k_chimera", "outcome")
    )

  S <- length(sample_names)

  summary_estimates <-
    outcome_counts %>%
    mutate(
      n_cugs = 1,
      n_chimeras = (k_chimera != 1)
    ) %>%
    summarize_at(
      vars(
        r,
        k_chimera,
        sample_names,
        n_chimeras,
        n_cugs
      ),
      list(~ sum(n * .))
    ) %>%
    rename(
      n_reads = "r",
      n_molecules = "k_chimera"
    ) %>%
    mutate(
      p_chimeras = n_chimeras / n_cugs,
      u = (n_molecules - n_cugs) / n_cugs,
      n_chimeras = NULL,
      n_molecules = NULL
    ) %>%
    select(n_cugs, p_chimeras, u, n_reads, everything())


  marginal <-
    summary_estimates %>%
    select(sample_names) %>%
    gather(sample, n_reads) %>%
    mutate(prop_reads = n_reads / sum(n_reads))


  marginal <-
    outcome_counts %>%
    summarize_at(
      sample_names,
      list(~ sum(n * as.integer(. > 0)))
    ) %>%
    data.table::transpose() %>%
    set_names("n_molecs") %>%
    bind_cols(marginal, .) %>%
    mutate(
      p_molecs = n_molecs / sum(n_molecs),
      RMR = n_reads / n_molecs
    )

  summary_estimates <-
    summary_estimates %>%
    select(n_cugs, p_chimeras, u, n_reads)

  conditional <-
    outcome_counts %>%
    ungroup() %>%
    mutate(
      n_reads = r,
      n_cugs = 1
    ) %>%
    select(r, n, n_reads, n_cugs, everything(), -outcome, -k_chimera) %>%
    group_by(r) %>%
    summarize_at(
      vars(n_reads, n_cugs, sample_names),
      list(~ sum(n * .))
    ) %>%
    mutate_at(
      vars(sample_names),
      list(~ . / n_reads)
    ) %>%
    ungroup() %>%
    mutate(
      p_cugs = n_cugs / sum(n_cugs),
      max_hop_rate = pmap_dbl(select(., sample_names), min) * (S - 1)
    ) %>%
    select(r, n_cugs, p_cugs, everything(), -n_reads)

  g <- get_prop_fugues(conditional$p_cugs, phat)


  summary_estimates <-
    summary_estimates %>%
    mutate(
      g = g,
      n_pm = round(n_cugs * (u + g)),
      RMR = n_reads / n_cugs
    ) %>%
    select(n_reads, n_cugs, n_pm, u, g, RMR, p_chimeras, everything())

  pi_r_hat <-
    estimate_pi_r_hat_matrix(
      conditional,
      sample_names,
      phat
    )

  return(
    list(
      summary_estimates = summary_estimates,
      marginal = marginal,
      conditional = conditional,
      pi_r_hat = pi_r_hat
    )
  )
}


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


create_outcome_counts <- function(read_counts, sample_names) {
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
      fill = list(nn = 0)
    ) %>%
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

fit_glm <- function(chimera_counts, S, max_r, conf_level = 0.99) {
  if (is.null(max_r)) {
    max_r <-
      chimera_counts %>%
      pull(r) %>%
      max()
  }


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
    unnest(c(
      tidied,
      confint_tidied,
      max_r
    )) %>%
    select(-std.error, -statistic, -p.value) %>%
    mutate_if(is.double,
      .funs = list(p = ~ exp(.))
    ) %>%
    select(
      max_r,
      estimate_p,
      conf.low_p,
      conf.high_p
    ) %>%
    rename(
      phat = estimate_p,
      phat_low = conf.low_p,
      phat_high = conf.high_p
    ) %>%
    mutate(
      SIHR = 1 - phat,
      SBIHR = (1 - phat) * ((4 * S - 1) / (4 * S - 4))
    )

  return(glm_estimates)
}


update_chimera_counts <- function(chimera_counts, glm_estimates) {
  chimera_counts <-
    chimera_counts %>%
    mutate(
      phat_chimeras = 1 - glm_estimates$phat^r,
      phat_chimeras_low = 1 - glm_estimates$phat_low^r,
      phat_chimeras_high = 1 - glm_estimates$phat_high^r
    )
  return(chimera_counts)
}


#' Step 4: estimate the sample index hopping rate (SIHR)
#' @param out A list *out* containing read_counts and sample_names.
#' @param max_r Maximum PCR duplication level to consider. The default considers all observed levels.
#' @return A list *out* with a total of 6 elements: read_counts, sample_names, outcome_counts, glm estimates, chimera counts, and summary_stats.
#' @details Given that the number of observations
#'  in the joined read counts datatable tends to be in the hundreds of millions,
#'  observations with similar outcomes are tallied to produce an outcome counts datatable which in turn
#'  is used to create a chimera counts datatable that tallies the number of chimeric and tho non-chimeric at each PCR duplication level r.
#'  The chimera counts datatable is then used as input to the GLM Model in order to estimate, from the proportion of chimeric observations,
#'  the sample index hopping rate. A list of summary statistics is also returned containing
#'  the conditional and marginal distributions of reads, and the
#'  proportion of molecules across samples \eqn{\pi_r} for each each PCR duplication level r.
#' @order 4
#' @export
estimate_hopping_rate <- function(out, max_r = NULL) {
  tic("Estimating SIHR. Step 1: creating outcome counts datatable")
  outcome_counts <-
    create_outcome_counts(out$read_counts, out$sample_names)
  toc()

  tic("Estimating SIHR. Step 2: creating chimera counts datatable")
  S <- length(out$sample_names)
  chimera_counts <- create_chimera_counts(outcome_counts, S)
  toc()

  tic("Estimating SIHR. Step 3: fitting GLM")
  glm_estimates <-
    fit_glm(chimera_counts,
      S = S,
      max_r
    )
  toc()

  tic("Estimating SIHR. Step 4: computing summary statistics")
  chimera_counts <-
    update_chimera_counts(
      chimera_counts,
      glm_estimates
    )

  summary_stats <- compute_summary_stats(
    outcome_counts,
    glm_estimates$phat
  )
  toc()

  out <- c(
    out,
    list(
      outcome_counts = outcome_counts,
      glm_estimates = glm_estimates,
      chimera_counts = chimera_counts,
      summary_stats = summary_stats
    )
  )
  return(out)
}
