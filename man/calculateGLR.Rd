\name{calculateGLR}
\alias{calculateGLR}
\title{Profile likelihood calculation using regression models}
\description{
  This is the function that calculates profileLikelihood for a single SNP. The main function \code{\link{evian}} calls this function repeatedly to obtain results for multiple SNPs.
}
\usage{
  calculateGLR(snp, formula_tofit, model, data, bim, lolim, hilim, m, bse,
                     family, c, plinkCC)
}
\arguments{
  \item{snp}{a string specifying the SNP of interests to be calculated.}
  \item{formula_tofit}{a formula object of the genetic model. The model should be formatted as \code{y~nuisance parameters}. The parameter of interest should not be included here.}
  \item{model}{a string specifying the mode of inheritance parameterization: \cr \code{additive, dominant, recessive, or overdominance}. See details.}
  \item{data}{data frame; read from the argument \code{data} in the main function \code{\link{evian}}. It should contain the SNP ID specified in the \code{snp} argument as a column name.}
  \item{bim}{data frame; read from from the argument \code{bim} in the main function \code{\link{evian}}. Provides allele information (base pair, effect/reference alleles) for the SNP of interest.}
  \item{lolim}{numeric; the lower limit for the grid or the minimum value of the regression parameter {\eqn{\beta}} used to calculate the likelihood function. }
  \item{hilim}{numeric; the upper limit for the grid or the maximum value of the regression parameter {\eqn{\beta}} used to calculate the likelihood funciton. }
  \item{m}{numeric; the density of the grid at which to compute the standardized likelihood function. A beta grid is defined as the grid of values for the SNP parameter used to evaluate the likelihood function.}
  \item{bse}{numeric; the number of beta standard errors to utilize in constraining the beta grid limits. Beta grid is evaluated at {\eqn{\beta}} +/- \code{bse}*s.e.}
  \item{family}{the link function for \code{glm}.}
  \item{c}{numeric; interval of Null Hypothesis to be tested.}
  \item{plinkCC}{A boolean type that specifies how case/control are coded. case/control were coded 1/0 if it is FALSE, and were coded 2/1 if TRUE.}
}
\details{
  \code{calculateGLR} conducts a likelihood ratio test for testing the SNP of interest. It uses the same numerical approach as the main function \link{calculateEvianMLE} to construct the likelihood function and it is then testing whether the effect of the SNP falls in an interval (-c, c) instead of testing whether the effect is 0 as in the \link{calculateEvianMLE}.
}
\value{
  This function outputs a dataframe that contains the summary statistics of the profile likelihood calculation. It contains the following columns:
      \itemize{
        \item \code{GLR}: the estimated generalized Likelihood ratio, a value smaller than 1 indicating in favor of the null hypothesis whereas a value greater than 1 indicating in favor of the alternative hypothesis.
        \item \code{boundary}: the boundary where null hypothesis is defined. i.e. the value c in (-c, c)
        \item \code{AF}: allele frequency for the effective allele
        \item \code{SNP}: SNP ID
        \item \code{bp}: base pair position from the \code{bim} input
        \item \code{effect}, \code{ref}: the effective allele and the other allele from the \code{bim} input
      }
  }

\references{
Bickel, D. R. (2012). ``The strength of statistical evidence for composite hypotheses: Inference to the best explanation.'' Statistica Sinica, 22, 1147-1198.

Zhang, Z., \& Zhang, B. (2013). ``A likelihood paradigm forclinical trials. Journal of Statistical Theory and Practice'', 7, 157-177.
}
\author{
  Dr. Lisa J Strug \email{lisa.strug@utoronto.ca}
}
\note{
  When \code{lolim} or \code{hilim} are NOT defined, then the boundaries of the beta grid will be determined by the default \code{bse=5}, or a \code{bse} defined by the user. Otherwise, the user can define the exact beta grid boundaries using \code{lolim} and \code{hilim}.

  In some cases the beta grid (using \code{bse} or \code{lolim},\code{hilim}) may need to be increased substantially (\code{bse} as large as 15) if covariates are present in the formula. This is automatically dealt by the current function, but contributes to longer computation time to find the appropriate ranges. Estimation may become inaccurate with large number of correlated covariates, which is a known limitation of profile likelihoods.
}


