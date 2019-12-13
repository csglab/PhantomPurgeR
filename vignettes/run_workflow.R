

#output_dir <- file.path("/home/rfarouni/Documents", "hiseq4000")
#joined_counts_filepath <-  file.path(output_dir, sprintf("%s_read_counts.rds",  dataset_name))



# #dir.create(output_dir)


#library(rhdf5)
#library(DropletUtils) # install but not load
library(tidyverse)
library(matrixStats)
library(broom)
library(furrr)
library(tictoc)
library(data.table)
library(cowplot)
library(scales)
library(zeallot) # For %<-% that unpacks lists in the Python manner
options(future.fork.enable =TRUE)

dataset_name <-"hiseq4000"
tic("Step 7: create umi counts matrices")

toc()


Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
```{Rcpp}
//[[Rcpp::depends(beachmat)]]
#include <beachmat/numeric_matrix.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace Rcpp;

// Checking for scalar inputs.

template<typename T, class V>
  T check_scalar(Rcpp::RObject incoming, const char* arg, const char* val) {
    V vec(incoming);
    if (vec.size()!=1) {
      std::stringstream err;
      err << arg << " should be " << val;
      throw std::runtime_error(err.str());
    }
    return vec[0];
  }

int check_integer_scalar(Rcpp::RObject incoming, const char* arg) {
  return check_scalar<int, Rcpp::IntegerVector>(incoming, arg, "an integer scalar");
}

double check_numeric_scalar(Rcpp::RObject incoming, const char* arg) {
  return check_scalar<double, Rcpp::NumericVector>(incoming, arg, "a numeric scalar");
}

bool check_logical_scalar(Rcpp::RObject incoming, const char* arg) {
  return check_scalar<bool, Rcpp::LogicalVector>(incoming, arg, "a logical scalar");
}



template<class V>
  std::vector<V> process_list(Rcpp::List incoming) {
    const size_t nsamples=incoming.size();
    std::vector<V> output(nsamples);
    for (size_t i=0; i<output.size(); ++i) {
      output[i] = V(incoming[i]);
    }
    return(output);
  }

template<class U, class V>
  void compare_lists(U left, V right) {
    if (left.size()!=right.size()) {
      throw std::runtime_error("lists are not of the same length");
    }
    const size_t nsamples=left.size();
    for (size_t i=0; i<nsamples; ++i) {
      if (left[i].size()!=right[i].size()) {
        throw std::runtime_error("list vectors are not of the same length");
      }
    }
    return;
  }

struct molecule {
  molecule (int s, size_t i, int g, int u) : index(i), sample(s), gene(g), umi(u) {}

  // need to handle situations where one sample has >2e9 UMIs.
  // otherwise, using ints for memory efficiency.
  size_t index;
  int sample, gene, umi;
};

/* Identifies which molecules should be retained in which samples,
* given the cell, gene and UMI combination for each molecule per sample.
* Also returns a diagnostic matrix of molecule-sample read counts.
*/
  //[[Rcpp::export]]
SEXP find_swapped2(SEXP cells, SEXP genes, SEXP umis, SEXP reads, SEXP minfrac, SEXP diagnostics) {
  BEGIN_RCPP

  auto Cells=process_list<Rcpp::StringVector>(cells);
  auto Genes=process_list<Rcpp::IntegerVector>(genes);
  auto Umis=process_list<Rcpp::IntegerVector>(umis);
  auto Reads=process_list<Rcpp::IntegerVector>(reads);

  compare_lists(Cells, Genes);
  compare_lists(Cells, Umis);
  compare_lists(Cells, Reads);

  const double mf=check_numeric_scalar(minfrac, "minimum fraction");
  const int diagcode=check_numeric_scalar(diagnostics, "diagcode");

  // Setting up the ordering vector.
  const size_t nsamples=Cells.size();
  size_t nmolecules=0;
  for (size_t i=0; i<nsamples; ++i) {
    nmolecules+=Cells[i].size();
  }

  std::vector<molecule> ordering;
  ordering.reserve(nmolecules);
  for (size_t i=0; i<nsamples; ++i) {
    const size_t cur_nmol=Cells[i].size();
    const auto& cur_genes=Genes[i];
    const auto& cur_umis=Umis[i];

    auto gIt=cur_genes.begin();
    auto uIt=cur_umis.begin();
    for (size_t j=0; j<cur_nmol; ++j, ++gIt, ++uIt) {
      ordering.push_back(molecule(i, j, *gIt, *uIt));
    }
  }

  // Sorting the indices.
  std::sort(ordering.begin(), ordering.end(), [&](const molecule& left, const molecule& right) {
    if (left.gene < right.gene) {
      return true;
    } else if (left.gene > right.gene) {
      return false;
    }

    if (left.umi < right.umi) {
      return true;
    } else if (left.umi > right.umi) {
      return false;
    }

    // Referencing the StringVectors to avoid copying them.
    return Cells[left.sample][left.index] < Cells[right.sample][right.index];
  });

  // Setting up the output, indicating which values to keep from each sample.
  std::vector<Rcpp::LogicalVector> notswapped(nsamples);
  for (size_t i=0; i<nsamples; ++i) {
    notswapped[i]=Rcpp::LogicalVector(Cells[i].size());
  }

  auto same_combination = [&] (const molecule& left, const molecule& right) {
    return left.gene==right.gene && left.umi==right.umi
    && Cells[left.sample][left.index]==Cells[right.sample][right.index];
  };

  // Iterating across runs of the same UMI/gene/cell combination.
  auto ostart=ordering.begin(), oend=ordering.begin();
  size_t nunique=0;

  while (ostart!=ordering.end()) {
    int max_nread=Reads[ostart->sample][ostart->index];
    int total_nreads=max_nread;
    auto best_mol=ostart;

    // ostart is always equal to oend at this point, so incrementing the latter to get to the next read.
    ++oend;
    while (oend!=ordering.end() && same_combination(*ostart, *oend)) {
      const int current_nread=Reads[oend->sample][oend->index];
      if (current_nread > max_nread) {
        max_nread=current_nread;
        best_mol=oend;
      }

      total_nreads += current_nread;
      ++oend;
    }

    if (double(max_nread)/total_nreads >= mf) {
      notswapped[best_mol->sample][best_mol->index]=1;
    }
    if (diagcode) {
      ++nunique;
    }
    ostart=oend;
  }

  // Creating the output list.
  Rcpp::List outlist(nsamples);
  for (size_t i=0; i<nsamples; ++i) {
    outlist[i]=notswapped[i];
  }
  Rcpp::List output(2);
  output[0]=outlist;
  output[1]=R_NilValue;

  // Storing diagnostic information about each unique combination.
  if (diagcode) {
    Rcpp::IntegerVector indices(nsamples);
    Rcpp::NumericVector values(nsamples);
    auto diag_out=beachmat::create_numeric_output(nunique, nsamples,
                                                  diagcode==1 ? beachmat::output_param("dgCMatrix", "Matrix") :
                                                    beachmat::output_param("HDF5Matrix", "HDF5Array"));

    auto ostart=ordering.begin(), oend=ordering.begin();
    size_t counter=0;
    while (ostart!=ordering.end()) {
      size_t nnzero=0;

      // Adding read counts per molecule, storing them in the matrix, then wiping them.
      while (oend!=ordering.end() && (ostart==oend || same_combination(*ostart, *oend))) {
        if (nnzero >= nsamples) {
          throw std::runtime_error("multiple instances of the same combination observed in a single sample");
        }
        indices[nnzero]=oend->sample;
        values[nnzero]=Reads[oend->sample][oend->index];
        ++oend;
        ++nnzero;
      }

      diag_out->set_row_indexed(counter, nnzero, indices.begin(), values.begin());
      ostart=oend;
      ++counter;
    }

    output[1]=diag_out->yield();
  }

  return output;
  END_RCPP
}
```
```{r}
sum(is.na(genes[[1]]))
sort(table(), decreasing = T)
map(genes, max)
```

```{r}
out$read_counts
```

```{r}

discard_nonmapped_counts <- function(cells, umis, genes, nreads)
  cur.cells <- cells[[i]]
cur.genes <- genes[[i]]
indx_keep <- cur.genes < length(unique(cur.genes))
cur.cells <- cur.cells[indx_keep]
cur.genes <- cur.genes[indx_keep])
```
cells = cells, umis = umis, genes = genes, nreads = nreads, ref.genes = ref.genes

```{r}
.check_indices <- function(val, tab, err) {
  if (is.character(val)) {
    if (is.null(tab)) {
      tab <- sort(unique(val))
    }
    val <- match(val, tab)
    if (any(is.na(val))) {
      stop(sprintf("entry of '%s' not in 'all.%ss'", err, err))
    }
    nvals <- length(tab)
  } else {
    if (!is.null(tab)) {
      nvals <- length(tab)
      if (length(val) && length(tab)<max(val)) {
        stop(sprintf("length of 'all.%ss' is less than 'max(%s)'", err, err))
      }
    } else {
      if (length(val)) {
        nvals <- max(val)
      } else {
        nvals <- 0L
      }
    }
    if (length(val) && min(val)<=0) {
      stop(sprintf("indices in '%s' must be positive", err))
    }
  }
  return(list(val=val, tab=tab, N=nvals))
}
```

```{r}
nsamples <- length(cells)
cleaned <- swapped <- vector("list", nsamples)
names(cleaned) <- names(swapped) <- names(cells)
i=1
for (i in seq_len(nsamples)) {



  all.cells <- sort(unique(cur.cells))

  stopifnot(identical(length(cur.cells), length(cur.genes)))



  if (is.null(value)) {
    value <- rep(1, length(cur.genes))
  }
  else {
    stopifnot(identical(length(cur.cells), length(value)))
  }
  gene.out <- .check_indices(cur.genes, all.genes, "gene")
  gene <- gene.out$val
  all.genes <- gene.out$tab
  ngenes <- gene.out$N
  cell.out <- .check_indices(cur.cells, all.cells, "cell")
  cell <- cell.out$val
  all.cells <- cell.out$tab
  ncells <- cell.out$N


  cleaned[[i]] <- makeCountMatrix(cur.genes[notswap], cur.cells[notswap],
                                  all.genes = ref.genes, all.cells = all.cells)
  if (get.swapped) {
    curswap <- !notswap
    swapped[[i]] <- makeCountMatrix(cur.genes[curswap],
                                    cur.cells[curswap], all.genes = ref.genes, all.cells = all.cells)
  }
}
```



```{r}
sum(rowSums(x) ==0)
```

```{r}
swap.out <- find_swapped2(cells, genes, umis, nreads, min.frac, 1)
```
```{r}
umi_counts_cell_gene$A1
```


```{r}
x<-as.matrix(swap.out[[2]])
```
```{r}
read_counts <- readRDS(joined_counts_filepath)
```

#' Join read counts of all samples
#' @param read_counts A list of read counts of all the samples
#' @param joined_counts_filepath If provided, the filepath of the joined readcounts
#' @return A joined dataframe of read counts and an outcome variable
#' @export
get_joined_read_counts <- function(samples, joined_counts_filepath = NULL) {
  if (!is.null(joined_counts_filepath) & file.exists(joined_counts_filepath)) {
    read_counts <- readRDS(joined_counts_filepath)
  } else if ((!is.null(joined_counts_filepath) & !file.exists(joined_counts_filepath)) | is.null(joined_counts_filepath)) {

  }


  return(read_counts)
}
if (!is.null(joined_counts_filepath)) {
  saveRDS(read_counts, joined_counts_filepath)
}
read_counts$gene <- ref_genes[read_counts$gene]



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

#' Load CellRanger 10xv3 data
#' @param sample The name of the filename (i.e. sample)
#' @return A dataframe
#' @importFrom rhdf5 h5read
read10xMolInfo_3 <- function(sample) {
  data <- list()
  all.barcodes <- as.vector(h5read(sample, "/barcodes"))
  all.barcodes <-
    sub("-[0-9]+", "", all.barcodes) # remove GEM group.
  data$cell <-
    all.barcodes[as.vector(h5read(sample, "/barcode_idx")) + 1L]
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/feature_idx")) + 1L
  data$reads <- as.vector(h5read(sample, "/count"))
  gene_name <- h5read(sample, "/features/name")
  gene_ids <- h5read(sample, "/features/id")


  return(data = data)
}

#' Load CellRanger 10xv2 data
#' The function uses cxx_get_cell_barcodes from dropUtils
#' @param sample The name of the filename (i.e. sample)
#' @return A dataframe
#' @importFrom rhdf5 h5read
read10xMolInfo_2 <- function(sample) {
  data <- list()

  data$cell <-
    .Call(
      DropletUtils:::cxx_get_cell_barcodes,
      sample,
      "barcode",
      NULL
    )
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/gene")) + 1L
  data$reads <- as.vector(h5read(sample, "/reads"))
  gene_name <- h5read(sample, "/gene_names")
  gene_ids <- h5read(sample, "/gene_ids")

  nonunique <-
    duplicated(gene_name) |
    duplicated(gene_name, fromLast = TRUE)

  gene_name[nonunique] <-
    str_c(gene_name[nonunique],
          gene_ids[nonunique],
          sep = ":"
    )

  data <- do.call(tibble, data)

  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene = gene_name[gene])

  return(data)
}

#' Load CellRanger data  for all samples
#' @param input_dir The filepath of the directory
#' @return A dataframe of read counts for all samples
#' @importFrom rhdf5 h5ls
#' @export
load_molecule_info_data <- function(input_dir) {
  molecule_info_filepaths <- get_h5_filenames(input_dir)

  available <-
    h5ls(molecule_info_filepaths[[1]],
         recursive = FALSE
    )
  version <- if ("barcode_idx" %in% available$name) {
    "3"
  } else {
    "2"
  }

  if (version == "3") {
    read_counts <- future_map(
      molecule_info_filepaths,
      read10xMolInfo_3
    )
  }
  else {
    read_counts <- map(
      molecule_info_filepaths,
      read10xMolInfo_2
    )
  }

  return(read_counts)
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

