% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob.dropout.R
\name{prob.dropout}
\alias{prob.dropout}
\title{Estimate drop-out probability of scRNA-seq expression data}
\usage{
prob.dropout(input, is.count = T, offsets = rep(0, ncol(input)),
  mcore = 1, ...)
}
\arguments{
\item{input}{matrix containing scRNA-seq expression values,
raw count or TPM values. Row is gene, col is sample.}

\item{is.count}{logical; if TRUE the input is raw count,
otherwise, it is TPM.}

\item{offsets}{offsets required for Poisson-NB model, the defaul is 
0 for all samples. To account for the influence of 
library size on the estimates of drop-out probability, 
user could use log transformed library size (total counts) 
as offsets.}

\item{mcore}{number of cores used in calculations. The parallel 
computation depends on package foreach and doParallel.}

\item{...}{further arguments passed to param.est and drop.prob functions.}
}
\value{
a matrix as the same dimension of the input and 
        containing drop-out probability for each gene at each sample
}
\description{
Estimate drop-out probability of scRNA-seq expression 
             data (raw count or TPM values) via mixture model. 
             Poisson-NB model for raw count; Gamma-Normal model for 
             TPM (or other normalized expression values).
}
\examples{
# raw count example
data(zebrafish)

offsets <- as.numeric(log(colSums(zebrafish)))
count <- zebrafish[rowSums(zebrafish > 5) > 4, ]

dp.mat <- prob.dropout(input = count, is.count = T, offsets = offsets, mcore = 3)


# TPM example
data(lung)

tpm <- lung[rowSums(lung>5) >4, ]

dp.mat <- prob.dropout(input = tpm, is.count = F, mcore = 3)

}
