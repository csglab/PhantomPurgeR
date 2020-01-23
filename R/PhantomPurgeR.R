#' PhantomPurgeR Workflow
#'
#' @details
#' See workflow vignettes on package's \href{https://csglab.github.io/PhantomPurgeR/pages/rpackage.html}{website}.
#'
#' Description of the workflow steps.
#' \enumerate{
#' \item  \code{\link{get_h5_filenames}}
#' \item  \code{\link{read10xMolInfoSamples}}
#' \item  \code{\link{join_data}}
#' \item  \code{\link{estimate_hopping_rate}}
#' \item  \code{\link{purge_phantoms}}
#' \item  \code{\link{make_plots}}
#' }
#' @docType package
#' @name PhantomPurgeR
#' @importFrom utils data
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

if (getRversion() >= "2.15.1")
  utils::globalVariables(c("."))
