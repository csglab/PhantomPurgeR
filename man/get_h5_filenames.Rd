% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_load_data.R
\name{get_h5_filenames}
\alias{get_h5_filenames}
\title{Step 1: return named list of sample filepaths}
\usage{
get_h5_filenames(input_dir)
}
\arguments{
\item{input_dir}{The directory in which the \emph{molecule_info.h5} files are located.}
}
\value{
A named list of samples filepaths.
}
\description{
Step 1: return named list of sample filepaths
}
\details{
The function read \emph{molecule_info.h5} files produced
by 10X Genomics CellRanger software found in user-provided \emph{input_dir}. Either rename the files in the \emph{input_dir} folder
as \emph{{sample_name}.h5} before or rename the named list after.
}
