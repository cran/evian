\name{evian}
\alias{evian}
\title{Evidential analysis for genetic data using regression models}
\description{
  Calculates the likelihood intervals for genetic association in a genomic region of interest. Covariates can be accommodated.
}
\usage{
evian(data, bim, xcols = NULL, ycol = NULL, covariateCol = NULL,
      formula = NULL, robust = FALSE, model='additive', m=200,
      bse = 5, lolim = NULL, hilim = NULL, kcutoff = c(8,32,100,1000),
      multiThread = 1, family='gaussian',plinkCC=F)
}
\arguments{
  \item{data}{a data frame includes a column for the response variable, one or multiple columns of genotype data (coded as \code{0, 1, 2,} or \code{NA}), and optionally columns for covariates. Headers are assumed. If the data is from related individuals, an additional column named \code{'FID'} needs to be included to specify the related structure. Using the PLINK toolkit with option \code{--recodeA} can produce the file in the required format and is recommended.}
 \item{bim}{a data frame with six columns representing chromosome, SNP ID, physical distance, base pair position, effective allele, and reference allele. i.e. data from a file in PLINK binary format (bim). No header is assumed, but the ordering of the columns must follow the standard bim file format.}
  \item{ycol}{numeric; column index in the \code{data} data frame for the column representing the response variable.}
  \item{xcols}{numeric vector; the column range in the \code{data} where genotype information is stored. Note that although a range of X is required, only one SNP at a time is calculated.}
  \item{covariateCol}{numeric or numeric vector; optional argument specifying which columns represent covariates. If left as \code{NULL}, no covariates will be included and the model \code{Y~snp} will be used.}
  \item{formula}{string; this is an alternative way of specifying model rather than using \code{xcols} and \code{ycol} arguments. This model follows the same format as the \code{glm} function (e.g. \code{Y~snp1+age+sex}). Note that in the case where multiple SNPs are included, only one SNP will be considered (e.g. given \code{Y~snp1+snp2}, the function will consider snp1 as the parameter of interests). The function can automatically identify SNPs with rsID as proper Xs, and would treat all other predictors as covariates.}
  \item{robust}{logical; default \code{FALSE}. If \code{TRUE}, then a robust adjustment is applied to the likelihood function to account for clustering in the data; See \link{robust_forCluster}.}
  \item{model}{a string that specifies the mode of inheritance parameterization: \cr \code{additive, dominant, recessive,} or \code{overdominance}. Default \code{additive}.}
  \item{m}{numeric; the density of the grid at which to compute the standardized likelihood function. A beta grid is defined as the grid of values for the SNP parameter used to
    evaluate the likelihood function.}
  \item{bse}{numeric; the number of beta standard errors to utilize in constraining the beta grid limits. Beta grid is evaluated at {\eqn{\beta}} +/- \code{bse}*s.e.}
  \item{lolim}{numeric; the lower limit for the grid or the minimum value of the regression parameter \eqn{\beta} used to calculate the likelihood function. }
  \item{hilim}{numeric; the upper limit for the grid or the maximum value of the regression parameter \eqn{\beta} used to calculate the likelihood funciton. }
  \item{kcutoff}{numeric or numeric vector; default = \code{c(8,32,100,1000)}. The strength of evidence criterion k. The function will calculate the \code{1/k} standardized likelihood intervals for each value provided.}
  \item{multiThread}{numeric; number of threads to use for parallel computing.}
  \item{family}{the link function for \code{glm}.}
  \item{plinkCC}{A boolean type that specifies how case/control are coded. case/control were coded 1/0 if it is FALSE, and were coded 2/1 if TRUE.}

}
\details{
 \code{evian} is the main function called to calculate the \code{1/k} likelihood intervals for the additive, dominant, recessive, or overdominance genotypic models. This function calls \code{calculateEvianMLE} in parallel to calculate the likelihood for each SNP. The calculation details can be found in \code{calculateEvianMLE}.

The input for the \code{data} and \code{bim} arguments can be obtained from the PLINK files; \code{data} is expected to follow PLINK format when run with the \code{--recodeA} option and \code{bim} can be obtained directly from a PLINK binary format file. Note if covariates are to be included, it is expected that the covariates are appended to the \code{data} file with a header for each covariate.

The statistical model can be specified in two ways. Column index can be provided through the \code{xcols}, \code{ycol}, and \code{covariateCol} arguments or through the \code{formula} argument, which can accept a formula specified as the \code{formula} argument in the R \code{glm} function. We recommend using \code{xcols}, \code{ycol}, and \code{covariateCol} arguments in most scenarios as this is relatively easier to input and it works for all the cases that we have considered so far. The alternative \code{formula} argument is not able to detect non-rsID variants as parameters of interests, and is only suggested in the scenario where only a single variant is of interest and that its rsID is known in advance. Since the profileLikelihood can only accomendate scalar parameter and thus if multiple rsID variants are inputted through \code{formula} option, it will only assume the first one to be parameter of interests.s

Parallel computing is avaliable through the use of the \code{multiThread} argument. This parallelization uses the \code{foreach} and \code{doMC} packages and will typically reduce computation time significantly. Due to this dependency, parallelization is not available on Windows OS as \code{foreach} and \code{doMC} are not supported on Windows.
}
\value{
This function outputs the row-combined the results from \code{calculateEvianMLE} for each of the SNPs included in the data/bim files. The exact output for each SNP can be found in the \code{calculateEvianMLE} documentation.
}

\note{
  When \code{lolim/hilim} is NOT defined, then the boundaries of the beta grid will be determined by the default \code{bse=5}, or a \code{bse} defined by the user. Otherwise, the user can define the exact beta grid boundaries using \code{lolim/hilim}.

  In some cases the beta grid (using \code{bse} or \code{lolim/hilim}) may need to be increased substantially (\code{bse} as large as 15) if covariates are present in the formula. This is automatically dealt by the current function, but contribute to longer computation time to find the appropriate ranges. Estimation may become inaccurate with large number of correlated covariates, which is a known limitation of profile likelihoods.
}
\seealso{
 \code{\link{calculateEvianMLE}}
}
\examples{

data(evian_linear_raw)
data(evian_linear_bim)
rst1=evian(data=evian_linear_raw, bim=evian_linear_bim, xcols=10:ncol(evian_linear_raw),
ycol=6, covariateCol=c(5,7:9), robust=FALSE, model="additive", m=200, lolim=-0.4,
hilim=0.4, kcutoff = c(32,100), multiThread=1,family='gaussian',plinkCC=FALSE)


}


