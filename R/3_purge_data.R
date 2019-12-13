
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


dedup_reads <- function(read_counts, sample_names) {
  read_counts <-
    read_counts %>%
    mutate_at(
      c(sample_names, paste0(sample_names, "_hat")),
      list(~ as.integer(. > 0))
    )

  return(read_counts)
}

#' make sparse count matrices
#' @param umi_counts list of sample data
#' @param all_genes ref_genes
#' @param return_discarded return discarded data
#' @return list of two lists containing retained and discarded data
#' @importFrom DropletUtils makeCountMatrix
make_count_matrices <- function(umi_counts, all_genes, return_discarded=TRUE ) {

  nsamples <- length(umi_counts)
  retained <- discarded <- vector("list", nsamples)
  names(retained) <- names(discarded) <- names(umi_counts)

  for (i in seq_len(nsamples)) {
    cur_sample <- umi_counts[[i]]

    cur_sample <-
      cur_sample %>%
      filter(retained +purged> 0 )

    cur_cells <- cur_sample$cell
    all_cells <- sort(unique(cur_cells))
    cur_genes <- cur_sample$gene
    cur_values_retained <- cur_sample$retained

    retained[[i]] <- makeCountMatrix(cur_genes,
                                     cur_cells,
                                     all.genes = all_genes,
                                     all.cells = all_cells,
                                     value=cur_values_retained)

    if (return_discarded) {
      cur_values_purged   <- cur_sample$purged
      discarded[[i]] <- makeCountMatrix(cur_genes,
                                        cur_cells,
                                        all.genes = all_genes,
                                        all.cells = all_cells,
                                        value=cur_values_purged)
    }

  }
  out <- list(retained = retained)
  if (return_discarded) {
    out$discarded <- discarded
  }
  return(out)
}


get_purge_summary <- function(umi_counts) {

  purge_summary <-
    umi_counts %>%
    select(c("sample", "purged_phantom", "purged_real", "retained")) %>%
    group_by(sample) %>%
    summarize_all(~sum(.) ) %>%
    mutate(umi_total= purged_phantom + purged_real + retained,
           phantom_prop= purged_phantom/umi_total) %>%
    select(sample, umi_total, retained, purged_real, purged_phantom, phantom_prop)
  return(purge_summary)

}

#' Purge and save read_counts
#' @param read_counts read counts
#' @param sample_names sample names
#' @return umi_counts dataset
create_umi_counts <- function(read_counts,
                              sample_names) {

  S <- length(sample_names)

  read_counts <- dedup_reads(read_counts, sample_names)

  read_counts <-
    read_counts %>%
    rename_at(sample_names, list(~ paste0("m", 1:S))) %>%
    rename_at(paste0(sample_names, "_hat"), list(~ paste0("t", 1:S)))


  # tallying molecule counts by cell-barcode and gene ID
  setDT(read_counts)
  umi_counts <-
    read_counts[, lapply(.SD, sum),
                keyby = "cell,gene,retain",
                .SDcols = -c("umi")
                ]

  # tranforming cell-gene molecule tally table into long format

  umi_counts <-
    melt(
      umi_counts,
      id = 1:3,
      variable.name = "sample",
      measure = patterns(m = "^m", t = "^t"),
      variable.factor = FALSE
    )

  sample_key <-
    tibble(
      sample = as.character(1:S),
      sample_name = sample_names
    )

  umi_counts <-
    left_join(umi_counts,
              sample_key,
              by = "sample"
    )

  umi_counts <-
    umi_counts %>%
    mutate(
      sample = sample_name,
      purged_phantom = m - t, # discarded phantom
      purged_real = t - t * retain, # discarded real
      retained = t * retain, # retained
      t = NULL,
      m = NULL,
      sample_name = NULL,
      retain = NULL
    )

  return(umi_counts )
}

split_counts_into_list <- function(umi_counts){
  sample_index <-
    umi_counts %>%
    pull(sample)

  umi_counts <-
    umi_counts %>%
    mutate(
      purged=purged_phantom+purged_real,
      purged_phantom=NULL,
      purged_real=NULL,
      sample = NULL
    )

  umi_counts <- split(
    umi_counts,
    sample_index
  )
  return(umi_counts)
}

merge_counts_outcomes <- function(read_counts, outcome_counts, sample_names){

  read_counts <-
    left_join(read_counts %>%
                select(outcome, cell, umi, gene, sample_names),
              outcome_counts %>%
                select(c("outcome", "retain", paste0(sample_names, "_hat"))),
              by = c("outcome")
    ) %>%
    select(-outcome)


  return(read_counts)
}

#' Reassign reads
#' @param out out from previous steps
#' @param torc TOR cutoff
#' @param return_readcounts If true the joined readcounts is returned
#' @param return_discarded return discarded data
#' @return out list with additional objects; umi_counts and summary_stats
#' @export
purge_phantoms <- function(out, torc, return_readcounts=FALSE, return_discarded=TRUE) {

  outcome_counts <- out$outcome_counts

  summary_stats <- compute_summary_stats(outcome_counts,
                                         out$glm_estimates$phat)

  outcome_counts <-
    reassign_reads(
      outcome_counts,
      summary_stats$pi_r_hat,
      out$glm_estimates$phat
    )

  outcome_counts <-
    mark_retained_observations(
      outcome_counts,
      summary_stats,
      torc
    )

  summary_stats <- get_threshold(outcome_counts, summary_stats)

  sample_names <- out$sample_names


  out$read_counts <- merge_counts_outcomes(out$read_counts, outcome_counts, sample_names)

  outcome_counts <-
    outcome_counts %>%
    arrange(-qr) %>%
    select(-c(paste0(sample_names, "_hat")))


  umi_counts <- create_umi_counts(out$read_counts, sample_names)

  summary_stats$purge_summary  <- get_purge_summary(umi_counts)

  umi_counts <- split_counts_into_list(umi_counts)

  message("making matrices")

  all_genes <- out$ref_genes

  umi_counts <- make_count_matrices(umi_counts, all_genes, return_discarded=return_discarded )

  out$outcome_counts <- outcome_counts

  if (!(return_readcounts)) {
    out$read_counts <- NULL
  }

  return(c( list(umi_counts=umi_counts, summary_stats=summary_stats), out))
}

#' Compute summary statistics
#' @param samples samples
#' @param barcode.length the length of the cell barcode (for v2 data)
#' @param max_r Maximum PCR duplication level to consider
#' @return list
#' @export
load_and_estimate_sihr <- function(samples,  barcode_length = NULL, max_r = NULL){

  out <- read10xMolInfoSamples(samples, barcode_length)

  out <- join_read_counts(out)

  out <- estimate_hopping_rate(out, max_r = NULL)

  return(out)
}


