
\name{glr_plot}
\alias{glr_plot}
\title{Plot methods for multiple glr result in a genomic region.}
\description{
  This function plots the log10 GLR value in the region calculated from the \code{\link{glr}} functions.
}
\usage{
 glr_plot(glr, snpName,kcut,col,legend_cex=1,...)
}
\arguments{
  \item{glr}{a dataframe, output from \code{\link{glr}}.}
  \item{snpName}{a string specifying the SNP to be marked.}
  \item{kcut}{numeric; the cut-off to be plotted. Can have value <=1.}
  \item{col}{color of the likelihood interval lines from smallest to largest. For instance, \code{c('red','green')} for LIs of \code{k=c(8,32)} means that the 1/8 interval will be plotted as red, and 1/32 will be plotted as green.}
  \item{legend_cex}{numeric; control the size of the legend, default 1.}
  \item{...}{Argument passed to the \code{plot} function.}
}
\details{
  This function takes output from \code{\link{glr}} as input. It will plot the GLR value on a log-scale and marked the SNPs of interests (defaultly the SNP with the maximum GLR values).
}


\examples{

  data(evian_linear_raw)
  data(evian_linear_bim)
  rst2=glr(data=evian_linear_raw, bim=evian_linear_bim, xcols=10:ncol(evian_linear_raw),
             ycol=6, covariateCol=c(5,7:9), c=0.025,robust=F, model="additive",
             m=200, lolim=-0.6, hilim=0.6, multiThread=1,family='gaussian',plinkCC=FALSE)

  # Plot the density for rs912
  glr_plot(glr=rst2,ylim=c(-0.5,2))


}


