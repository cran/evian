% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{eviandata_linear}
\alias{eviandata_linear}
\title{Example dataset with a quantative outcome.}
\description{
  This dataset included the genotypic and phenotypic information of 1444 individuals in proper formats. This dataset is used together with \code{\link{evianmap_linear}} to illustrate how to use the main function \code{\link{evian_linear}}.
}
\usage{
data(eviandata_linear)
}

\details{
  This is an example dataset for \code{\link{evian_linear}} function. It contained 1444 individuals, and for each of the individuals, their genotype at 10 SNPs and a continuous outcome (Y_norma) were stored. Three additional covariates (Fev, BMI_group, Age_group) were provided as well. Specifically, our function can incorporate with related individuals, some of these individuals in the dataset are correlated with others, and are specified through the \code{FID} column.
}
