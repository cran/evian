% --- Source file: ./man/evian_linear.Rd ---
% Part of the evian package
% Distributed under GPL 2 or later

\name{evian_logit}
\alias{evian_logit}
\title{Evidential analysis for binary outcome data using logistic regression models}
\description{
  Calculates the likelihood intervals for genetic association of a binary trait in a genomic region of interest. Covariates can be accommodated.
}
\usage{
evian_logit(data, bim, xcols = NULL, ycol = NULL, covariateCol = NULL, formula = NULL,
  robust = FALSE, model = "additive", m = 200, bse = NULL, lolim = log(0.025),
  hilim = log(4), kcutoff = c(8,32,100,1000), multiThread = 1)
}
\arguments{
  \item{data}{a data frame includes a column for the response variable, one or multiple columns of genotype data (coded as \code{0, 1, 2,} or \code{NA}), and optionally columns for covariates. Headers are assumed. If the data is from related individuals, an additional column named \code{'FID'} needs to be included to specify the related structure. Using the PLINK toolkit with option \code{--recodeA} can produce the file in the required format and is recommended. Note specifically for binary outcome, we assumed that the case/control status were coded as \code{0/1} rather than \code{1/2} as how PLINK codes them.}
 \item{bim}{data frame; six columns representing chromosome, SNP ID, physical distance, base pair position, effective allele, and reference allele. i.e. data from a file in PLINK binary format (bim). No header is assumed, but the order must follow the standard bim file format.}
  \item{ycol}{numeric; column index in the \code{data} data frame for the column representing the response variable.}
  \item{xcols}{numeric vector; the column range in the \code{data} where genotype information is stored. Note that although a range of X is required, only one SNP at a time is calculated.}
  \item{covariateCol}{numeric or numeric vector; optional argument specifying which columns represent covariates. If left as \code{NULL}, no covariates will be included and the model \code{Y~snp} will be used.}
  \item{formula}{string; this is an alternative way of specifying model rather than using \code{xcols} and \code{ycol} arguments. This model follows the same formate as the \code{glm} function (e.g. \code{Y~snp1+age+sex}). Note that in the case where multiple SNPs are included, only one SNP will be calculated at a time (e.g. given \code{Y~snp1+snp2}, the function will estimate \code{Y~snp1} and \code{Y~snp2} separately). The function can automatically identify snps with rsID as proper Xs, and would treat all other predictors as covariates.}
  \item{robust}{logical; default \code{FALSE}. If \code{TRUE}, then a robust adjustment is applied to the likelihood function to account for clustering in the data; See \link{robust_forCluster}.}
  \item{model}{string; specifies the mode of inheritance parameterization: \cr \code{additive, dominant, recessive,} or \code{overdominance}. Default \code{additive}.}
  \item{m}{numeric; the density of the grid at which to compute the standardized likelihood function. A beta grid is defined as the grid of values for the SNP parameter used to
    evaluate the likelihood function.}
  \item{bse}{numeric; the number of beta standard errors to utilize in constraining the beta grid limits. Beta grid is evaluated at {\eqn{\beta}} +/- \code{bse}*s.e. The default for this option is NULL, meaning to use \code{lolim} and \code{hilim} values specified. See details \code{\link{calculateLogitMLE}}.}
  \item{lolim}{numeric; the lower limit for the grid or the minimum value of the regression parameter \eqn{\beta} used to calculate the likelihood function. }
  \item{hilim}{numeric; the upper limit for the grid or the maximum value of the regression parameter \eqn{\beta} used to calculate the likelihood funciton. }
  \item{kcutoff}{numeric or numeric vector; default = \code{c(8,32,100,1000)}. The strength of evidence criterion k. The function will calculate the \code{1/k} standardized likelihood intervals for each value provided.}
  \item{multiThread}{numeric; number of threads to use for parallel computing.}
}
\details{
 \code{evian_logistic} is the main function called to calculate the \code{1/k} likelihood intervals for the additive, dominant, recessive, or overdominance genotypic models when a binary phenotype is presented. This function calls \code{calculateLogitMLE} in parallel to calculate the likelihood for each SNP. The calculation details can be found in \code{calculateLogitMLE}.

The inputs for the \code{data} and \code{bim} arguments can be obtained from the PLINK files; \code{data} is expected to follow PLINK format when run with the \code{--recodeA} option and \code{bim} can be obtained directly from a PLINK binary format file. Note if covariates are to be included, it is expected that the covariates are appended to the \code{data} file with a header for each covariate.

The statistical model can be specified in two ways. Column index can be provided through the \code{xcols}, \code{ycol}, and \code{covariateCol} arguments or through the \code{formula} argument, which can accept a formula specified as the \code{formula} argument in the R \code{glm} function. We recommend using \code{xcols}, \code{ycol}, and \code{covariateCol} arguments in most scenarios as this is relatively easier to input and it works for all the cases that we have considered so far. The alternative \code{formula} argument is not able to detect non-rsID variants as parameters of interests, and is only suggested in the scenario where only a few variants are of interest and all of their rsID are known in advance. In this case, using this argument can save the time to search for the column index values corresponding to those SNPs.

Parallel computing is avaliable through the use of the \code{multiThread} argument. This parallelization uses the \code{foreach} and \code{doMC} packages and will typically reduce computation time significantly. Due to this dependency, parallelization is not available on Windows OS as \code{foreach} and \code{doMC} are not supported on Windows.
}
\value{
This function outputs the row-combined the results from \code{calculateLogitMLE} for each of the SNPs included in the data/bim files. The exact output for each SNP can be found in the \code{calculateLogitMLE} documentation.
}

\note{
  In some cases the beta grid (using \code{bse} or \code{lolim/hilim}) may need to be increased substantially (\code{bse} as large as 15) if covariates are present in the formula. This is automatically dealt by the current function, but contribute to longer computation time to find the appropriate ranges. Estimation may become inaccurate with large number of correlated covariates, which is a known limitation of profile likelihoods.
}
\seealso{
 \code{\link{calculateLogitMLE}}
}
\examples{

data(eviandata_binary)
data(evianmap_binary)
\dontshow{
rst1=evian_logit(data=eviandata_binary, bim=evianmap_binary, xcols=10:14, ycol=6, robust=FALSE,
model="additive", m=100, kcutoff = c(32,100), multiThread=1)
}
\donttest{
rst1=evian_logit(data=eviandata_binary, bim=evianmap_binary, xcols=10:19, ycol=6, robust=FALSE,
model="additive", m=1000, kcutoff = c(32,100), multiThread=1)

#Alternatively you can use the formula argument to run the same model as above

rst2=evian_logit(data=eviandata_binary, bim=evianmap_binary, formula='PHENOTYPE~rs461+rs462+rs477+
rs479+rs491+rs492+rs504+rs509+rs519+rs542', robust=FALSE, model="additive",
m=1000, kcutoff = c(32,100), multiThread=1)
}
}



