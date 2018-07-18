% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{eviandata_binary}
\alias{eviandata_binary}
\title{Example dataset with a binary outcome.}
\description{
  This dataset included the genotypic and phenotypic information of 250 individuals in proper formats. This dataset is used together with \code{\link{evianmap_binary}} to illustrate how to use the main function \code{\link{evian_logit}}.
}
\usage{
data(eviandata_binary)
}

\details{
  This is an example dataset for \code{\link{evian_logit}} function. It contained 250 individuals, and for each of the individuals, their genotype at 30 SNPs and a binary outcome (PHENOTYPE; coded as 0/1) were stored. Three additional covariates (age, weight, city) were provided as well. Specifically, our function can incorporate with related individuals, some of these individuals in the dataset are correlated with others, and are specified through the \code{FID} column.
}
