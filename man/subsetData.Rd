% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{subsetData}
\alias{subsetData}
\title{interior subsetting function}
\description{
  This is an internal function that subsets the SNP column with matching name and removes rows with missing observations.
}
\usage{
  subsetData(snp,formula_tofit,data,plinkCC)
}
\arguments{
  \item{snp}{a string specifying the SNP of interests. The SNP ID must exist in \code{data}.}
  \item{formula_tofit}{a formula object of the genetic model. This is directly obtained from \code{\link{evian}} function.}
  \item{data}{a data frame inputted from the main function. Should contain the SNP ID \code{snp} as one of the column names.}
  \item{plinkCC}{A boolean type that specifies how case/control are coded. case/control were coded 1/0 if it is FALSE, and were coded 2/1 if TRUE.}
}
\details{
  \code{subsetData} is an interior function that subsets the full dataset into a smaller set containing only one specific SNP by the \code{snp} option. It will then remove any rows with missing values.
}
\value{
This function returned a dataframe containing phenotype, covariates in their original column names as in the full dataset, and a column called \code{X} representing the genotype information for the SNP chosen. The column names are essetial.
}

\keyword{models}

