% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wrap.R
\name{Impute}
\alias{Impute}
\title{Dropout imputation using gene-specific best-performing methods}
\usage{
\code{Impute(data, method.choice, scale = 1, pseudo.count = 1,
count_path, labels = NULL, cell.clusters = NULL, drop_thre = NULL, cores = 4,
cluster.type = "SOCK", type = "TPM", network.path = "../network/network.rds",
...)}
}
\arguments{
\item{data}{matrix; raw counts (genes as rows and samples as columns)}

\item{method.choice}{character; best performing method in training data for
each gene}

\item{scale}{integer; scaling factor to divide all expression levels by
(defaults to 1)}

\item{pseudo.count}{integer; pseudo-count to be added to expression levels
to avoid log(0) (defaults to 1)}

\item{count_path}{character; path to data file}

\item{labels}{character; vector specifying the cell type of each column of
\code{data}}

\item{cell.clusters}{integer; number of cell subpopulations}

\item{drop_thre}{numeric; between 0 and 1 specifying the threshold to
determine dropout values}

\item{cores}{integer; number of cores used for paralell computation}

\item{cluster.type}{character; either "SOCK" or "MPI"}

\item{network.path}{character; path to .rds file containing network
coefficients and intercept}
}
\value{
matrix; imputed and normalized expression values
}
\description{
\code{Impute} performs dropout imputation based on the
performance results obtained in the training data, coupled to normalization
using \code{normalization.function}
}
\details{
Values that are 0 in \code{data} are imputed according to the
best-performing methods indicated in \code{methods}.
}
\seealso{
\code{\link{EvaluateMethods}},
\code{\link{ImputeBaseline}},
\code{\link{ImputeNetwork}},
\code{\link{ImputeScImpute}},
\code{\link{scimpute}}
}
