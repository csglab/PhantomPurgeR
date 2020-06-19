
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
    arrange(-qr, tor) %>%
    slice(1:2) %>%
    arrange(-tor) %>%
    mutate(approach = c("torc_before", "discard_torc"))


  next_thresh <-
    outcome_counts %>%
    filter(!retain) %>%
    arrange(qr, tor) %>%
    slice(1) %>%
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
      FN = round(get_product_logsum(pmax(1 - o + FPp - g, 1e-16), n_cugs)),
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


make_count_matrices <- function(umi_counts, all_genes, return_discarded = TRUE) {
  nsamples <- length(umi_counts)
  retained <- discarded <- vector("list", nsamples)
  names(retained) <- names(discarded) <- names(umi_counts)

  for (i in seq_len(nsamples)) {
    cur_sample <- umi_counts[[i]]

    cur_sample <-
      cur_sample %>%
      filter(retained + purged > 0)

    cur_cells <- cur_sample$cell
    all_cells <- sort(unique(cur_cells))
    cur_genes <- cur_sample$gene
    cur_values_retained <- cur_sample$retained

    retained[[i]] <- makeCountMatrix(cur_genes,
      cur_cells,
      all.genes = all_genes,
      all.cells = all_cells,
      value = cur_values_retained
    )

    if (return_discarded) {
      cur_values_purged <- cur_sample$purged
      discarded[[i]] <- makeCountMatrix(cur_genes,
        cur_cells,
        all.genes = all_genes,
        all.cells = all_cells,
        value = cur_values_purged
      )
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
    summarize_all(~ sum(.)) %>%
    mutate(
      umi_total = purged_phantom + purged_real + retained,
      phantom_prop = purged_phantom / umi_total
    ) %>%
    select(sample, umi_total, retained, purged_real, purged_phantom, phantom_prop)
  return(purge_summary)
}


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

  return(umi_counts)
}

split_counts_into_list <- function(umi_counts) {
  sample_index <-
    umi_counts %>%
    pull(sample)

  umi_counts <-
    umi_counts %>%
    mutate(
      purged = purged_phantom + purged_real,
      purged_phantom = NULL,
      purged_real = NULL,
      sample = NULL
    )

  umi_counts <- split(
    umi_counts,
    sample_index
  )
  return(umi_counts)
}

merge_counts_outcomes <- function(read_counts, outcome_counts, sample_names) {
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

#' Step 5: Reassign hopped reads and purges phantom molecules
#' @param out out from previous steps
#' @param torc TOR cutoff
#' @param return_readcounts If true the joined readcounts is returned
#' @param return_discarded return discarded data
#' @return A the initial list *out* with umi_counts, summary_stats added.
#' @details For each outcome, the conditional posterior probability \eqn{q|y} of the possible true
#'  samples of origin is computed  by plugging in  \eqn{\pi_r}, the estimated
#'  proportion of molecules across samples. The index of the sample with the maximum posterior
#'  probability along with posterior probability itself is added to the original joined read count table.
#'  The predicted true sample of origin and its associated posterior probability is then used to reassign
#'  reads to their predicted sample of origin.
#'
#'  In order to remove predicted phantom molecules from the data while minimizing the rate of false positives and false negatives,
#'  the Trade-Off Ratio (TOR) statistic is computed by dividing the marginal increase in FNs over the marginal decrease
#'  of FPs for each observed unique \eqn{qr*} value. The cutoff \eqn{TOR*} that gets effectively chosen would correspond to the largest
#'  observed TOR value not exceeding the preset TOR cutoff value (default value is 3). All molecules with corresponding TOR values
#'  strictly less than  \eqn{TOR*} cutoff - not TOR cutoff- are discarded. For example,  if we have \eqn{tor= (0.1, 0.5, 2.9, 4.1, ...)}
#'  and TOR cutoff=3, then \eqn{TOR*} cutoff=2.9 and predicted real molecules corresponding to \eqn{tor=0.1} and \eqn{tor=0.5} are discarded.
#'  To purge the data, the read counts are first deduplicated to obtain a table of molecule (i.e. UMI) counts.
#'
#'  After purging, the molecule counts are collapsed over gene labels to produce a gene-by-cell umi-count expression matrices
#'  for all the samples sequenced in the same lane.
#' @order 5
#' @export
purge_phantoms <- function(out, torc, return_readcounts = FALSE, return_discarded = TRUE) {
  sample_names <- out$sample_names
  outcome_counts <- out$outcome_counts
  summary_stats <- out$summary_stats

  tic("Purging phantom molecules. Step 1: reassigning hopped reads")

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
  toc()

  tic("Purging phantom molecules. Step 2: getting observed tor threshold below user-provided cutoff")

  summary_stats <- get_threshold(outcome_counts, summary_stats)
  toc()

  tic("Purging phantom molecules. Step 3: marking retained observations in read counts table")

  out$read_counts <- merge_counts_outcomes(
    out$read_counts,
    outcome_counts,
    sample_names
  )
  toc()


  outcome_counts <-
    outcome_counts %>%
    arrange(-qr) %>%
    select(-c(paste0(sample_names, "_hat")))


  umi_counts <- create_umi_counts(out$read_counts, sample_names)

  summary_stats$purge_summary <- get_purge_summary(umi_counts)



  tic("Purging phantom molecules. Step 4: creating sparse count matrices of cleaned and discarded data")

  umi_counts <- split_counts_into_list(umi_counts)


  all_genes <- out$ref_genes

  umi_counts <- make_count_matrices(umi_counts,
    all_genes,
    return_discarded = return_discarded
  )

  out$summary_stats <- summary_stats
  out$outcome_counts <- outcome_counts

  toc()

  if (!(return_readcounts)) {
    out$read_counts <- NULL
  }

  return(c(list(umi_counts = umi_counts), out))
}



#' A function wraper for running the PhantomPurgeR workflow in a single step.
#' @param samples A named list of filepaths
#' @param torc TOR cutoff
#' @param max_r Maximum PCR duplication level to consider
#' @param barcode_length the length of the cell barcode for v2 data
#' @param min_umi_cell The minimum number of UMIs associated with a cell barcode. Barcodes with fewer UMIs are discarded.
#'  Default is 1 (i.e. all cell barcodes are retained).
#'  For large datasets with many samples set to a larger value (50-200) to reduce memory requirements.
#' @param return_readcounts If true the joined readcounts is returned
#' @param return_discarded return discarded data
#' @return out list pf umi_counts, metadata, and summary statistics
#' @export
phantom_purger <- function(samples, torc, max_r = NULL, barcode_length = NULL, min_umi_cell=1, return_readcounts = FALSE, return_discarded = TRUE) {
  out <- read10xMolInfoSamples(samples, barcode_length = barcode_length)
  out <- join_data(out, min_umi_cell)
  out <- estimate_hopping_rate(out, max_r = max_r)
  out <- purge_phantoms(out,
    torc = torc,
    return_readcounts = return_readcounts,
    return_discarded = return_discarded
  )

  return(out)
}
