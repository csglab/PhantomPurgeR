
plan(multiprocess)
dataset_name <-"hiseq4000"
input_dir <- file.path("/home/rfarouni/Documents/index_hopping/data/hiseq4000/input")
output_dir <- file.path("/home/rfarouni/Documents", dataset_name)
# #dir.create(output_dir)
joined_counts_filepath <-  file.path(output_dir, sprintf("%s_read_counts.rds",  dataset_name))
#
#read_counts <- load_molecule_info_data(input_dir)
# read_counts <- get_joined_counts_table(read_counts, joined_counts_filepath = joined_counts_filepath)
read_counts <- get_joined_counts_table( joined_counts_filepath = joined_counts_filepath)
read_counts <- add_outcome_variable(read_counts)
outcome_counts <- create_outcome_counts(read_counts)
fit_out <- estimate_hopping_rate(outcome_counts)

c(read_counts, outcome_counts, summary_stats) %<-%
  reassign_reads_and_mark_retained_observations(read_counts,
    outcome_counts,
    fit_out,
    torc=3
  )


toc()

tic("Step 6: Purge and save read counts datatable to disk")



purge_and_save_read_counts(
  read_counts,
  dataset_name,
  summary_stats$sample_names,
  output_dir
)

toc()


tic("Step 7: create umi counts matrices")
umi_counts_cell_gene <-
  create_umi_counts(
    read_counts,
    sample_names
  )
toc()



read_counts <-
  read_counts %>%
  select(-c("retain", paste0(sample_names, "_hat")))


data_list <-
  list(
    umi_counts_cell_gene = umi_counts_cell_gene,
    read_counts = read_counts,
    outcome_counts = outcome_counts,
    fit_out = fit_out,
    summary_stats = summary_stats
  )

toc()
