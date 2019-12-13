#' Workflow: Phantom Purge
#' \enumerate{
#' \item  \code{\link{load_molecule_info_data}}
#' loads molecule_info files and creates read counts datatable.
#' Each sample's *molecule_info.h5* file should be renamed to
#'  *{sample_name}.h5* and placed in the input_dir folder.
#' \item  \code{\link{get_joined_counts_table}}
#' join read counts of all samples
#' \item  \code{\link{create_outcome_counts}}
#' create outcome counts datatable
#' \item  \code{\link{estimate_hopping_rate}}
#' create a chimera counts datatable and estimate hopping rate
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
