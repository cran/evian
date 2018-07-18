% --- Source file: ./man/robust_for_cluster.Rd ---
% Part of the evian package
% Distributed under GPL 2 or later

\name{robust_forCluster}
\alias{robust_forCluster}
\title{Robust adjustment function}
\description{
The robust function computes an adjustment that is applied to the likelihood function to account for the cluster nature of the data.
}
\usage{
robust_forCluster(formula, data, family)
}
\arguments{
  \item{formula}{a formula specifying the response and possible covariates to keep in the output data frame. This is directly obtained from the \code{evian_linear} or \code{evian_logistic} function.}
  \item{data}{data frame; from the output of \code{subsetData}.}
  \item{family}{string; the link function for \code{glm}. Currently this only supported linear and logistic.}
}
\details{
The robust function is called from within \code{evian_logit} or \code{evian_linear} functions. It computes a robust adjustment factor that is
applied to the likelihood function to account for the cluster nature of the data. The family ID column (FID) specifies the clusters. The robust adjustment factor is the ratio of the regular variance estimator of the maximum likelihood estimate (MLE) to the sandwich variance estimator of the MLE, where the `meat' of the sandwich variance estimator is corrected for clustering in the data (Blume et.al, 2007).
If the data is not clustered (i.e. the observations are independent) then the adjustment factor can still be applied to make the working model robust to possible model misspecifications (Royall and Tsou, 2003).
}
\value{
A numeric constant.
}

\references{
Blume, J. D., Su, L., Remigio, M. O., & McGarvey, S. T. (2007). Statistical evidence for GLM regression parameters: A robust likelihood approach. Statistics in Medicine, 26, 2919-2936.

Royall, R. , Tsou, T. S. (2003). Interpreting statistical evidence by using imperfect models: robust adjusted likelihood functions. J Roy Stat Soc B; 65: 391-404.
}

\author{ Zeynep Baskurt \email{zeynep.baskurt@sickkids.ca} }


