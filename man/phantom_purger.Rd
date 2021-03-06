% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_purge_data.R
\name{phantom_purger}
\alias{phantom_purger}
\title{A function wraper for running the PhantomPurgeR workflow in a single step.}
\usage{
phantom_purger(
  samples,
  torc,
  max_r = NULL,
  barcode_length = NULL,
  min_umi_cell = 1,
  return_readcounts = FALSE,
  return_discarded = TRUE
)
}
\arguments{
\item{samples}{A named list of filepaths}

\item{torc}{TOR cutoff}

\item{max_r}{Maximum PCR duplication level to consider}

\item{barcode_length}{the length of the cell barcode for v2 data}

\item{min_umi_cell}{The minimum number of UMIs associated with a cell barcode. Barcodes with fewer UMIs are discarded.
Default is 1 (i.e. all cell barcodes are retained).
For large datasets with many samples set to a larger value (50-200) to reduce memory requirements.}

\item{return_readcounts}{If true the joined readcounts is returned}

\item{return_discarded}{return discarded data}
}
\value{
out list pf umi_counts, metadata, and summary statistics
}
\description{
A function wraper for running the PhantomPurgeR workflow in a single step.
}
