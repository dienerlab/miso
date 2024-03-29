% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_counts}
\alias{plot_counts}
\title{Plots counts for several taxa across a co-variable}
\usage{
plot_counts(
  ps,
  variable,
  tax_level = "genus",
  taxa = NULL,
  normalized = TRUE,
  pc = 0.5,
  only_data = FALSE,
  zeros = TRUE
)
}
\arguments{
\item{ps}{A phyloseq data set.}

\item{variable}{The name of the co-variable.}

\item{tax_level}{The taxonomy level to use. Defaults to genus.}

\item{taxa}{A character vector denoting the taxa to be plotted. Defaults
to plotting all taxa.}

\item{normalized}{Whether to normalize the counts using the DESeq2 size
factors.}

\item{pc}{The pseudo count to add.}

\item{only_data}{Only get the raw data for the plot.}

\item{zeros}{Whether to also include zero counts.}
}
\value{
Nothing.
}
\description{
Plots counts for several taxa across a co-variable
}
\examples{
 NULL

}
