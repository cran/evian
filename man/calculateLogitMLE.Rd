% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later


\name{calculateLogitMLE}
\alias{calculateLogitMLE}
\title{Profile likelihood calculation for binary outcome data using logistic regression models}
\description{
  This is the function that calculates profileLikelihood for a single SNP. The main function \code{\link{evian_logit}} calls this function repeatedly to obtain results for multiple SNPs.
}
\usage{
  calculateLogitMLE(snp,formula_tofit,model,data,bim,lolim,hilim,m,bse,k,robust,family)
}
\arguments{
  \item{snp}{a string specifyingthe SNP to be calculated.}
  \item{formula_tofit}{a formula object of the genetic model. The model should be formatted as \code{y~nuisance parameters}. The parameter of interest should not be included here.}
  \item{model}{a string specifying the mode of inheritance parameterization: \cr \code{additive, dominant, recessive, or overdominance}. See \code{\link{calculateLinearMLE}} for details.}
   \item{data}{data frame; read from the argument \code{data} in the main function \code{\link{evian_logit}}. It should contain the SNP ID specified in the \code{snp} argument as a column name.}
  \item{bim}{data frame; read from from the argument \code{bim} in the main function \code{\link{evian_logit}}. Provides allele information (base pair, effect/reference alleles) for the SNP of interest.}
  \item{lolim}{numeric; the lower limit for the grid or the minimum value of the regression parameter {\eqn{\beta}} used to calculate the likelihood function. }
  \item{hilim}{numeric; the upper limit for the grid or the maximum value of the regression parameter {\eqn{\beta}} used to calculate the likelihood funciton. }
  \item{m}{numeric; the density of the grid at which to compute the standardized likelihood function. A beta grid is defined as the grid of values for the SNP parameter used to
    evaluate the likelihood function.}
  \item{bse}{numeric; the number of beta standard errors to utilize in constraining the beta grid limits. Beta grid is evaluated at {\eqn{\beta}} +/- \code{bse}*s.e.}
  \item{k}{numeric or numeric vector; The strength of evidence criterion k. Reads from the input of \code{kcutoff} from the main function}
  \item{robust}{logical; if \code{TRUE}, then a robust adjustment is applied to the likelihood function to account for the cluster nature in the data. See \code{\link{robust_forCluster}}.}
  \item{family}{the link function for \code{glm}. Currently only linear (`gaussian') and logistic (`binomial') are supported. This is currently auto-filled by the main \code{\link{evian_logit}} function.}
}
\details{
  \code{calculateLogitMLE} calculates the profile likelihood for a single SNP. A proper grid range is first established for \eqn{\beta} then the standardized profile likelihood is evaluated at each of the \code{m} cuts uniformly spread across the grid. The choice of grid range can be either specified through \code{lolim} and \code{hilim} arguments (default) or can be specified through providing a numerical \code{bse} argument. This will enable the grid search feature as in \code{\link{evian_linear}}. After the standardized profile likelihood is evalualted, the MLE for \eqn{\beta} is computed as well as the likelihood intervals for each value of k provided.

Note when the \code{model} is \code{dominant} or \code{recessive}, the odds ratio cannot simply be inverted to obtain its reciprocal model, otherwise the heteozygous genotype group will be misclassified. Instead the genotype coding is swapped so that it codes with respect to the risk allele prior to running the analysis so odds ratio are properly calculated.}
\value{

  This function outputs a list containg 4 elements that can be directly accessed using \code{$} operator.
  \item{theta}{numeric vector; the \code{m} odds ratios used to estimate the standardized profile likelihood.}
  \item{profile.lik.norm}{numeric vector; the corresponding \code{m} standardized profile likelihood value at each of the odds ratios in \code{theta}. If \code{robust}=\code{TRUE}, then the values will be adjusted by the robust factor.}
  \item{k_cutoff}{numeric vector; It specifies which k-cutoff had been used in the calculation, ordered from the smallest k to the largest k.}
  \item{SummaryStats}{data frame; contains contains the summary statistics of the profile likelihood calculation:
      \itemize{
        \item \code{mle}: the estimates for SNP effect (odds Ratio for logistic, and \eqn{\beta} for linear) with respect to the effective allele
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
\author{
  Dr. Lisa J Strug \email{lisa.strug@utoronto.ca}
}
\note{
  As mentioned in the linear case, finding a proper beta grid range can be time-consuming. For logistic regression specifically, we assume that most of the odds ratios are between 0.025 to 4, and we assign those values as default for the lower and upper bound. If the grid search feature used in the linear case is preferred, specify a value for \code{bse} and leave \code{lolim} and \code{hilim} unchanged.
}


