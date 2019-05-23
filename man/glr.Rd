\name{glr}
\alias{glr}
\title{Generalized Likelihood Ratio test}
\description{
  Conducts a generalized likelihood ratio test testing whether {\eqn{\beta \in}  (-\code{c},\code{c})}. Covariates can be accommodated.
}
\usage{
  glr(data, bim, xcols = NULL, ycol = NULL, covariateCol = NULL,
      formula = NULL, c, robust = FALSE, model = 'additive',
      m = 200, bse = 5, lolim = NULL, hilim = NULL, multiThread = 1,
      family='gaussian',plinkCC=F)
}
\arguments{
  \item{data}{a data frame includes a column for the response variable, one or multiple columns of genotype data (coded as \code{0, 1, 2,} or \code{NA}), and optionally columns for covariates. Headers are assumed. If the data is from related individuals, an additional column named \code{'FID'} needs to be included to specify the related structure. Using the PLINK toolkit with option \code{--recodeA} can produce the file in the required format and is recommended.}
  \item{bim}{a data frame with six columns representing chromosome, SNP ID, physical distance, base pair position, effective allele, and reference allele. i.e. data from a file in PLINK binary format (bim). No header is assumed, but the ordering of the columns must follow the standard bim file format.}
  \item{ycol}{numeric; column index in the \code{data} data frame for the column representing the response variable.}
  \item{xcols}{numeric vector; the column range in the \code{data} where genotype information is stored. Note that although a range of X is required, only one SNP at a time is calculated.}
  \item{covariateCol}{numeric or numeric vector; optional argument specifying which columns represent covariates. If left as \code{NULL}, no covariates will be included and the model \code{Y~snp} will be used.}
  \item{formula}{string; this is an alternative way of specifying model rather than using \code{xcols} and \code{ycol} arguments. This model follows the same format as the \code{glm} function (e.g. \code{Y~snp1+age+sex}). Note that in the case where multiple SNPs are included, only the first SNP will be considered (e.g. given \code{Y~snp1+snp2}, the function will consider snp1 as the parameter of interests). The function can automatically identify SNPs with rsID as proper Xs, and would treat all other predictors as covariates.}
  \item{c}{numeric; interval of Null Hypothesis to be tested.}
  \item{robust}{logical; default \code{FALSE}. If \code{TRUE}, then a robust adjustment is applied to the likelihood function to account for clustering in the data; See \link{robust_forCluster}.}
  \item{model}{a string that specifies the mode of inheritance parameterization: \cr \code{additive, dominant, recessive,} or \code{overdominance}. Default \code{additive}.}
  \item{m}{numeric; the density of the grid at which to compute the standardized likelihood function. A beta grid is defined as the grid of values for the SNP parameter used to
    evaluate the likelihood function.}
  \item{bse}{numeric; the number of beta standard errors to utilize in constraining the beta grid limits. Beta grid is evaluated at {\eqn{\beta}} +/- \code{bse}*s.e.}
  \item{lolim}{numeric; the lower limit for the grid or the minimum value of the regression parameter \eqn{\beta} used to calculate the likelihood function. }
  \item{hilim}{numeric; the upper limit for the grid or the maximum value of the regression parameter \eqn{\beta} used to calculate the likelihood funciton. }
  \item{multiThread}{numeric; number of threads to use for parallel computing.}
  \item{family}{the link function for \code{glm}.}
  \item{plinkCC}{A boolean type that specifies how case/control are coded. case/control were coded 1/0 if it is FALSE, and were coded 2/1 if TRUE.}

}
\details{
  This is a similar function to the main function \code{evian}. Instead of testing H0: \eqn{\beta}=0, it tests for whether H0:{\eqn{\beta \in} (-\code{c},\code{c})}.
}
\value{
  This function outputs the row-combined the results from \code{calculateGLR} for each of the SNPs included in the data/bim files. The exact output for each SNP can be found in the \code{calculateGLR} documentation.
}

\note{
  When \code{lolim/hilim} is NOT defined, then the boundaries of the beta grid will be determined by the default \code{bse=5}, or a \code{bse} defined by the user. Otherwise, the user can define the exact beta grid boundaries using \code{lolim/hilim}.

  In some cases the beta grid (using \code{bse} or \code{lolim/hilim}) may need to be increased substantially (\code{bse} as large as 15) if covariates are present in the formula. This is automatically dealt by the current function, but contribute to longer computation time to find the appropriate ranges. Estimation may become inaccurate with large number of correlated covariates, which is a known limitation of profile likelihoods.
}
\seealso{
  \code{\link{calculateGLR}}
}
\examples{

  data(evian_linear_raw)
  data(evian_linear_bim)
  rst2=glr(data=evian_linear_raw, bim=evian_linear_bim, xcols=10:ncol(evian_linear_raw),
             ycol=6, covariateCol=c(5,7:9), c=0.025,robust=F, model="additive",
             m=200, lolim=-0.6, hilim=0.6, multiThread=1,family='gaussian',plinkCC=FALSE)


}


