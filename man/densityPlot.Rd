% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{densityPlot}
\alias{densityPlot}
\title{Plot profile likelihood density for a single SNP.}
\description{
  This function plots the density distribution for a single SNP calculated from the \code{\link{evian}} functions.
}
\usage{
  densityPlot(dList, snpName, kcut = NULL, pl = 'linear', xlim = NULL,
  color = c('red','orange','green','blue'), round = 2,legend_cex=1)
}
\arguments{
  \item{dList}{a row-combined list, output from \code{\link{evian}}.}
  \item{snpName}{a string specifying the SNP to be plotted.}
  \item{kcut}{numeric; the cut-off to be plotted. If \code{kcut}=\code{NULL}, all values of k in the \code{kcutoff} will be plotted.}
  \item{pl}{a string specifying the y-axis for the plot. The y-axis will be plotted as 'Odds Ratio' if \code{pl} is specified as \code{logit}, 'Beta' otherwise.}
  \item{xlim}{graphical parameter used in \code{plot} function}
  \item{color}{color of the likelihood interval lines from smallest to largest. For instance, \code{c('red','green')} for LIs of \code{k=c(8,32)} means that the 1/8 interval will be plotted as red, and 1/32 will be plotted as green.}
  \item{round}{numeric; number of digits displayed on the plot.}
  \item{legend_cex}{numeric; control the size of the legend, default 1.}
}
\details{
  This function takes output from \code{\link{evian}} as input. It will plot the density of the estimated standardized profile likelihood for the SNP of interest. Some basic summary statistics will be included on the plot too.
}

\references{
Strug, L. J., Hodge, S. E., Chiang, T., Pal, D. K., Corey, P. N., & Rohde, C. (2010). A pure likelihood approach to the analysis of genetic association data: an alternative to Bayesian and frequentist analysis. Eur J Hum Genet, 18(8), 933-941. doi:10.1038/ejhg.2010.47

Royall, R. (1997). Statistical Evidence: A Likelihood Paradigm. London, Chapman and Hall.
}

\examples{

data(evian_linear_raw)
data(evian_linear_bim)
rst1=evian(data=evian_linear_raw, bim=evian_linear_bim, xcols=10:ncol(evian_linear_raw),
      ycol=6, covariateCol=c(5,7:9), robust=FALSE, model="additive", m=200, lolim=-0.4,
      hilim=0.4, kcutoff = c(32,100), multiThread=1,family='gaussian',plinkCC=FALSE)

# Plot the density for rs912
densityPlot(dList=rst1,snpName='rs912')


}


