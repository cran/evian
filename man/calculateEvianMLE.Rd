\name{calculateEvianMLE}
\alias{calculateEvianMLE}
\title{Profile likelihood calculation using regression models}
\description{
  This is the function that calculates profileLikelihood for a single SNP. The main function \code{\link{evian}} calls this function repeatedly to obtain results for multiple SNPs.
}
\usage{
  calculateEvianMLE(snp, formula_tofit, model, data, bim, lolim, hilim,
                    m, bse, k, robust, family, plinkCC)
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
  \item{k}{numeric or numeric vector; The strength of evidence criterion k. Reads from the input of \code{kcutoff} from the main \code{\link{evian}} function}
  \item{robust}{logical; if \code{TRUE}, then a robust adjustment is applied to the likelihood function to account for the cluster nature in the data. See \link{robust_forCluster}.}
  \item{family}{the link function for \code{glm}.}
  \item{plinkCC}{A boolean type that specifies how case/control are coded. case/control were coded 1/0 if it is FALSE, and were coded 2/1 if TRUE.}
}
\details{
  \code{calculateEvianMLE} calculates the profile likelihood for a single SNP. A proper grid range is first established for \eqn{\beta} then the standardized profile likelihood is evaluated at each of the \code{m} cuts uniformly spread across the grid. Based on the standardized profile likelihood, the MLE for \eqn{\beta} is computed as well as the likelihood intervals for each value of \code{k} provided. \cr
  For different genetic models, their coding schemes are shown as below:
  \preformatted{
   Additive
 AA  0
 AB  1
 BB  2

   Dominant
 AA  0
 AB  1
 BB  1

   Recessive
 AA  0
 AB  0
 BB  1

 Overdominance model
        A  D
    AA  0  0
    AB  1  1
    BB  2  0
}
Specifically for the overdominance model, the column of interest is the D column.
}
\value{
  This function outputs a list containg 4 elements that can be directly accessed using '\code{$}' operator.
  \item{theta}{numeric vector; It stores all \code{m} \eqn{\beta} values that used to estimate the standardized profile likelihood.}
  \item{profile.lik.norm}{numeric vector; the corresponding \code{m} standardized profile likelihood value at each of the \eqn{\beta} values in \code{theta}. If \code{robust}=\code{TRUE}, then the values will be adjusted by the robust factor.}
  \item{k_cutoff}{numeric vector; It specifies which k-cutoff had been used in the calculation, ordered from the smallest k to the largest k.}
  \item{SummaryStats}{data frame; contains the summary statistics of the profile likelihood calculation. It contains the following columns:
      \itemize{
        \item \code{mle}: the estimates for SNP effect with respect to the effective allele
        \item \code{maxlr}: maximum likelihood ratio in the beta grid defined by \code{lolim} and \code{hilim}
        \item \code{AF}: allele frequency for the effective allele
        \item \code{SNP}: SNP ID
        \item \code{bp}: base pair position from the \code{bim} input
        \item \code{effect}, \code{ref}: the effective allele and the other allele from the \code{bim} input
        \item \code{robustFactor}: robust factor calculated, set to 1 if \code{robust}=\code{FALSE}.
        \item \code{lo_1}, \code{hi_1}, \code{lo_2}, \code{hi_2}...: the lower and upper bound of the likelihood intervals for the kth cut-off in \code{k_cutoff}.
      }
  }
}
\references{
 Strug, L. J., Hodge, S. E., Chiang, T., Pal, D. K., Corey, P. N., & Rohde, C. (2010). A pure likelihood approach to the analysis of genetic association data: an alternative to Bayesian and frequentist analysis. Eur J Hum Genet, 18(8), 933-941. doi:10.1038/ejhg.2010.47

 Strug, L. J., & Hodge, S. E. (2006). An alternative foundation for the planning and evaluation of linkage analysis. I. Decoupling "error probabilities" from "measures of evidence". Hum Hered, 61(3), 166-188. doi:10.1159/000094709

 Royall, R. (1997). Statistical Evidence: A Likelihood Paradigm. London, Chapman and Hall.
}

\note{
  When \code{lolim} or \code{hilim} are NOT defined, then the boundaries of the beta grid will be determined by the default \code{bse=5}, or a \code{bse} defined by the user. Otherwise, the user can define the exact beta grid boundaries using \code{lolim} and \code{hilim}.

  In some cases the beta grid (using \code{bse} or \code{lolim},\code{hilim}) may need to be increased substantially (\code{bse} as large as 15) if covariates are present in the formula. This is automatically dealt by the current function, but contributes to longer computation time to find the appropriate ranges. Estimation may become inaccurate with large number of correlated covariates, which is a known limitation of profile likelihoods.
}


