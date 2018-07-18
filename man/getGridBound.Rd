% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{getGridBound}
\alias{getGridBound}
\title{Obtain the range of the grid where MLE will be searched at}
\description{
  This is an internal function that provides the range where the profileLikelihood function would search for MLE.
}
\usage{
  getGridBound(formula,data,bse,k,m,family,robust)
}
\arguments{
  \item{formula}{a formula specifying the response and possible covariates to keep in the output dataframe. This is directly obtained from \code{\link{evian_linear}} or the \code{\link{evian_logit}} function.}
  \item{data}{a data frame inputted from the output of \code{\link{subsetData}}.}
  \item{bse}{numeric. The number of beta standard errors to utilize in constraining the beta grid limits. Passed down from argument \code{bse} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{k}{numeric vector. The strength of evidence criterion k. Passed down from argument \code{kcutoff} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{m}{numeric. The density of the grid at which to compute the standardized likelihood function. Passed down from argument \code{m} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{family}{a string representing the link function for \code{ProfileLikelihood::ProfileLikelihood.glm}. Currently only supports linear (`gaussian') and logistic (`binary').}
  \item{robust}{A numeric value, robust correction factor.}
}
\details{
  \code{getGridBound} is an interior function that searches for the proper grid range that would be used to search for MLE. This is done through two steps: First, it finds a starting grid range by fitting a (generalized) linear model to obtain the estimate and s.e. of the beta. Then the starting grid range can be defined as mean +- bse*s.e. . In the case where robust correction is needed, the grid will be defined as mean +- bse*s.e./correction factor. Then the function determines an optimal grid range by using \code{\link{expandBound}} function.
}
\value{
  This function returns a numeric vector of length 2 that represents the lower and upper bounds of the grid for the MLE search.
}

\keyword{models}

