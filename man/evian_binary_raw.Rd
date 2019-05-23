% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{evian_binary_raw}
\alias{evian_binary_raw}
\title{Example dataset with a binary outcome.}
\description{
  This dataset included the genotypic and phenotypic information of 250 individuals in proper formats. This dataset is used together with \code{\link{evian_binary_bim}} to illustrate how to use the main function \code{\link{evian}}.
}
\usage{
data(evian_binary_raw)
}

\details{
  This is an example dataset for \code{\link{evian}} function. It contained 250 individuals, and for each of the individuals, their genotype at 30 SNPs and a binary outcome (PHENOTYPE; coded as 0/1) were stored. Three additional covariates (age, weight, city) were provided as well. Specifically, our function can incorporate with related individuals, some of these individuals in the dataset are correlated with others, and are specified through the \code{FID} column. This is usually what a plink .raw file looks like except for the case/control status not coded as 1/2.
}
