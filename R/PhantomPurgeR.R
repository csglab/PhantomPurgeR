#' Workflow: Phantom Purge
#' \enumerate{
#' \item  \code{\link{get_h5_filenames}}
#' Returns named list of sample filepaths of CellRanger molecule_info.h5 files in provided directory.
#' \item  \code{\link{read10xMolInfoSamples}}
#' Loads molecule_info.h5 samples. Renamed files in the input_dir folder as
#'  *{sample_name}.h5* or rename the list of samples' filepaths.
#' \item  \code{\link{join_data}}
#' Join read counts of all samples. Merge mapped data by cell, umi, and gene keys.
#' \item  \code{\link{estimate_hopping_rate}}
#' Create outcome counts datatable, chimera counts datatable, and estimate the sample index hopping rate.
#' \item  \code{\link{purge_data}}
#' Reassign hopped reads and purge phantom molecules.
#' \item  \code{\link{make_plots}}
#' Create diagnostics plots.
#' }
#' @docType package
#' @name PhantomPurgeR
#' @import data.table
#' @import ggplot2
#' @import cowplot
#' @import tictoc
#' @importFrom matrixStats rowSums2 rowCounts
#' @importFrom DropletUtils read10xMolInfo makeCountMatrix
#' @importFrom broom confint_tidy tidy
#' @importFrom scales scientific
#' @importFrom tidyr nest unnest spread complete gather
#' @importFrom tibble enframe tibble
#' @importFrom purrr map2 map %>% set_names reduce pmap_dbl
#' @importFrom stringr str_c
#' @importFrom furrr future_map future_pmap_dfr
#' @importFrom dplyr left_join group_by ungroup mutate mutate_if mutate_at select filter summarize rename_at
#'  summarize_at summarize_all bind_cols bind_rows add_tally tally slice pull rename vars top_n n arrange everything matches
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
