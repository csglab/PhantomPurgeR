
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

#' Compute summary statistics
#' @param outcome_counts outcome dataset
#' @param max_r Maximum PCR duplication level to consider
#' @return list containing two dataframes (glm estimates and chimera counts)
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
      FRM = n_reads / n_molecs # ,
      #  n_molecs = n_molecs / 1000000,
      #   n_reads = n_reads / 1000000
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
      frm = n_reads / n_cugs
    ) %>%
    select(n_reads, n_cugs, n_pm, u, g, frm, p_chimeras, everything())

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

get_product_logsum <- function(x, y) {
  z_log <- log(x) + log(y)
  z <- exp(z_log)
  return(z)
}

infer_sample_of_origin_outcome <- function(..., log_zi, S, phat) {
  y <- c(...)[1:S]
  sample_pi <- c(...)[(S + 1):(2 * S)]
  log_sample_pi <- log(sample_pi)
  posterior_s <- y * log_zi + log_sample_pi

  ## normalize to the sample with maximum posterior probability
  posterior_s <- posterior_s - max(posterior_s)
  posterior_s <- exp(posterior_s)

  posterior <- posterior_s / sum(posterior_s)
  q <- max(posterior)
  s <- which.max(posterior)

  return(list(
    s = s,
    q = q
  ))
}


infer_sample_of_origin <- function(outcome_counts, pi_r_hat, S, phat) {
  outcome_counts <-
    left_join(
      outcome_counts,
      pi_r_hat,
      by = c("r"),
      suffix = c(".count", ".pi")
    )

  # zi <- ((S - 1) / (1 / p - 1))
  log_zi <- log((S - 1)) - log(1 / phat - 1)


  posteriors <-
    future_pmap_dfr(
      outcome_counts %>%
        select(matches(".count|.pi$")),
      log_zi = log_zi,
      S = S,
      phat = phat,
      infer_sample_of_origin_outcome
    )

  outcome_counts <-
    bind_cols(outcome_counts, posteriors) %>%
    select(-matches(".pi$|.count$"))

  return(outcome_counts = outcome_counts)
}



reassign_reads_outcome <- function(outcome_counts, S, sample_names) {
  outcome_counts_reassigned <-
    outcome_counts %>%
    mutate(
      rr = r,
      ss = s
    ) %>%
    spread(ss, rr, fill = 0L) %>%
    rename_at(vars(as.character(1:S)), list(~ paste0(sample_names, "_hat"))) %>%
    select(-k_chimera) %>%
    select(outcome, s, q, n, everything()) %>%
    arrange(q)

  return(outcome_counts_reassigned)
}




reassign_reads <- function(outcome_counts,
                           pi_r_hat,
                           phat) {

  sample_names <-
    setdiff(
      colnames(outcome_counts),
      c("n", "r", "k_chimera", "outcome")
    )
  S <- length(sample_names)

  outcome_counts <-
    infer_sample_of_origin(
      outcome_counts,
      pi_r_hat,
      S,
      phat
    )

  outcome_counts <- reassign_reads_outcome(outcome_counts, S, sample_names)

  return(outcome_counts)
}

get_threshold <- function(outcome_counts, summary_stats) {
  u <-
    summary_stats$summary_estimates %>%
    pull(u)

  g <-
    summary_stats$summary_estimates %>%
    pull(g)

  n_cugs <- summary_stats$summary_estimates$n_cugs

  current_thresh <-
    outcome_counts %>%
    filter(retain) %>%
    top_n(2, qr) %>%
    mutate(approach = c("torc_before", "discard_torc"))


  next_thresh <-
    outcome_counts %>%
    filter(!retain) %>%
    top_n(-1, qr) %>%
    mutate(approach = c("torc_after"))

  no_thresh <-
    outcome_counts %>%
    slice(n()) %>%
    mutate(approach = c("no_discarding"))


  cutoff_dt <-
    bind_rows(list(current_thresh, next_thresh, no_thresh)) %>%
    select(
      approach, outcome, s, qr, n, tor, o,
      FP, FN, TP, TN, FPm, FNm
    )

  no_purging <- list(
    "no_purging", NA, NA, NA, NA, NA, 1 + u,
    round(get_product_logsum(u + g, n_cugs)), 0,
    round(get_product_logsum(1 - g, n_cugs)), 0,
    NA, NA
  )

  names(no_purging) <- colnames(cutoff_dt)

  cutoff_dt <-
    bind_rows(cutoff_dt, no_purging)

  summary_stats$cutoff_dt <- cutoff_dt

  return(summary_stats)
}



mark_retained_observations <- function(outcome_counts, summary_stats, torc) {
  u <-
    summary_stats$summary_estimates %>%
    pull(u)

  g <-
    summary_stats$summary_estimates %>%
    pull(g)

  n_cugs <- summary_stats$summary_estimates$n_cugs


  outcome_counts <-
    outcome_counts %>%
    group_by(r) %>%
    add_tally(n, name = "nn") %>%
    ungroup() %>%
    mutate(
      p_outcome = n / sum(n)
    ) %>%
    select(-nn)

  outcome_counts <-
    outcome_counts %>%
    mutate(
      qr = 1 - q,
      qs = -10 * log10(1e-16) + 10 * log10(qr + 1e-16),
      q = NULL
    ) %>%
    arrange(qr) %>%
    mutate(
      o = cumsum(p_outcome),
      FPp = cumsum(get_product_logsum(
        qr,
        p_outcome
      )),
      FP = round(get_product_logsum(FPp, n_cugs)),
      FN = round(get_product_logsum(1 - o + FPp - g, n_cugs)),
      TP = round(get_product_logsum(o - FPp, n_cugs)),
      TN = round(get_product_logsum(u + g - FPp, n_cugs)),
      FPm = last(FP) - FP, # marginal
      FNm = FN - last(FN), # marginal
      tor = FNm / FPm
    )

  tor_thresh <-
    outcome_counts %>%
    summarize(tor_thresh = tor[which.max(-pmax(tor, torc))]) %>%
    pull(tor_thresh)

  outcome_counts <-
    outcome_counts %>%
    mutate(
      retain = (tor >= tor_thresh) & FPm > 1
    ) %>%
    select(
      outcome, s, qr, r, n, o,
      FP, FN, TP, TN,
      FPm, FNm, tor, retain, qs,
      FPp, p_outcome, everything()
    )

  return(outcome_counts)
}


#' Reassign reads
#' @param read_counts read counts
#' @param outcome_counts outcome dataset
#' @param fit_out Output of estimate_hopping_rate function
#' @param torc TOR cutoff
#' @return outcome dataset
#' @export
reassign_reads_and_mark_retained_observations <- function(read_counts, outcome_counts, fit_out, torc) {

  summary_stats <- compute_summary_stats(outcome_counts,
                                         fit_out$glm_estimates$phat)

  outcome_counts <-
    reassign_reads(
      outcome_counts,
      summary_stats$pi_r_hat,
      fit_out$glm_estimates$phat
    )

  outcome_counts <-
    mark_retained_observations(
      outcome_counts,
      summary_stats,
      torc
    )

  summary_stats <- get_threshold(outcome_counts, summary_stats)

  sample_names <-
    setdiff(
      colnames(outcome_counts),
      c("n", "r", "k_chimera", "outcome")
    )

  summary_stats$sample_names <- sample_names


  read_counts <-
    left_join(read_counts %>%
                select(outcome, cell, umi, gene, sample_names),
              outcome_counts %>%
                select(c("outcome", "retain", paste0(sample_names, "_hat"))),
              by = c("outcome")
    ) %>%
    select(-outcome)

  outcome_counts <-
    outcome_counts %>%
    arrange(-qr) %>%
    select(-c(paste0(sample_names, "_hat")))

  return(list(read_counts=read_counts, outcome_counts=outcome_counts, summary_stats=summary_stats))
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
