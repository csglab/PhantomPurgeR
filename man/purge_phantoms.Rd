% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_purge_data.R
\name{purge_phantoms}
\alias{purge_phantoms}
\title{Step 5: Reassign hopped reads and purges phantom molecules}
\usage{
purge_phantoms(out, torc, return_readcounts = FALSE, return_discarded = TRUE)
}
\arguments{
\item{out}{out from previous steps}

\item{torc}{TOR cutoff}

\item{return_readcounts}{If true the joined readcounts is returned}

\item{return_discarded}{return discarded data}
}
\value{
A the initial list \emph{out} with umi_counts, summary_stats added.
}
\description{
Step 5: Reassign hopped reads and purges phantom molecules
}
\details{
For each outcome, the conditional posterior probability \eqn{q|y} of the possible true
samples of origin is computed  by plugging in  \eqn{\pi_r}, the estimated
proportion of molecules across samples. The index of the sample with the maximum posterior
probability along with posterior probability itself is added to the original joined read count table.
The predicted true sample of origin and its associated posterior probability is then used to reassign
reads to their predicted sample of origin.

In order to remove predicted phantom molecules from the data while minimizing the rate of false positives and false negatives,
the Trade-Off Ratio (TOR) statistic is computed by dividing the marginal increase in FNs over the marginal decrease
of FPs for each observed unique \eqn{qr*} value. The cutoff \eqn{TOR*} that gets effectively chosen would correspond to the largest
observed TOR value not exceeding the preset TOR cutoff value (default value is 3). All molecules with corresponding TOR values
strictly less than  \eqn{TOR*} cutoff - not TOR cutoff- are discarded. For example,  if we have \eqn{tor= (0.1, 0.5, 2.9, 4.1, ...)}
and TOR cutoff=3, then \eqn{TOR*} cutoff=2.9 and predicted real molecules corresponding to \eqn{tor=0.1} and \eqn{tor=0.5} are discarded.
To purge the data, the read counts are first deduplicated to obtain a table of molecule (i.e. UMI) counts.

After purging, the molecule counts are collapsed over gene labels to produce a gene-by-cell umi-count expression matrices
for all the samples sequenced in the same lane.
}
