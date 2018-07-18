\name{expandBound}
\alias{expandBound}
\title{A recursive function that expands the grid search for MLE.}
\description{
  This is an internal function that finds the proper boundary of the grid.
}
\usage{
  expandBound(data,bse,parameters,formula,m,k,family)
}
\arguments{
  \item{data}{a data frame inputted from the main function.}
  \item{bse}{numeric. The number of beta standard errors to utilize in constraining the beta grid limits. Passed down from argument \code{bse} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{parameters}{a numeric vector of length 3 providing the starting values for the search. This is obtained from the \code{\link{getGridBound}} function. The three numeric values in the vector should represent the beta estimates, s.e., and the correction factor respectively. Details can be found in \code{\link{getGridBound}}.}
  \item{formula}{a formula specifying the response and possible covariates to keep in the output dataframe. This is directly obtained from \code{\link{evian_linear}} or the \code{\link{evian_logit}} function.}
  \item{k}{numeric vector. The strength of evidence criterion k. Passed down from argument \code{kcutoff} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{m}{numeric. The density of the grid at which to compute the standardized likelihood function. Passed down from argument \code{m} in the main \code{\link{evian_linear}} or \code{\link{evian_logit}} function.}
  \item{family}{a string representing the link function for \code{ProfileLikelihood::ProfileLikelihood.glm}. Currently only supports linear (`gaussian') and logistic (`binary').}
}
\details{
  Even though the initial grid bound calculated from \code{getGridBound} works for most of the data, there can be cases where \code{bse} needs to be increased in order to observe all the Likelihood Intervals (LIs) specified from the main function in the range \code{kcutoff} calculated. In this case, our approach is to check whether the current grid range includes the largest LIs. The function will expand the grid range by increasing \code{bse} by 1 if it is not included. This step will be running recursively until the largest LIs are included in the current grid range.
}
\value{
This function returns a numeric vector of length two representing the optimal lower and upper bounds for the grid on which the later functions will search for MLE.
}

\keyword{models}
