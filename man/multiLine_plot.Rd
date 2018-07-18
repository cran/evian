% --- Source file: ./man/evian_linear.Rd ---
  % Part of the evian package
% Distributed under GPL 2 or later

\name{multiLine_plot}
\alias{multiLine_plot}
\title{Plot methods for multiple likelihood intervals in a genomic region.}
\description{
  This function plots the likelihood intervals (LIs) for all SNPs calculated using either \code{evian_linear} or \code{evian_logit}.
}
\usage{
  multiLine_plot(bpstart = 0, bpend = 1000000000, dList, title = "My model",
  showmaxlr = 3, kcut = NULL, pl = 'linear', ylim = c(-0.5,10),
  color = c('violet','green','red','blue'), markSNP = NULL, round = 2)
}
\arguments{
  \item{bpstart,bpend}{numeric; indicating the range of base pairs to be plotted. From \code{bpstart} to \code{bpend}.}
  \item{dList}{a row-combined list, output from \code{evian_linear} or \code{evian_logit.}}
  \item{title}{string; title of plot}
  \item{showmaxlr}{numeric; number of top SNPs to display on the graph. Default = 3. SNPs are chosen by their maximum likelihood ratio values.}
  \item{kcut}{numeric; the cut-off to be plotted. If \code{kcut}=\code{NULL}, all intervals will be plotted.}
  \item{pl}{a string specifying the y-axis for the plot. The y-axis will be plotted as 'Odds Ratio' if \code{pl} is specified as \code{logit}, 'Beta' otherwise.}
  \item{markSNP}{vector of strings; indicates which SNPs to be marked on the plot. By default it will mark all SNPs that are significant at the smallest cut-off.}
  \item{round}{numeric; number of digits displayed on the plot.}
  \item{ylim}{graphical parameter used in \code{plot} function}
  \item{color}{color of the likelihood interval lines from smallest to largest. For instance, \code{c('red','green')} for LIs of \code{k=c(8,32)} means that the 1/8 interval will be plotted as red, and 1/32 will be plotted as green.}
}
\details{
  This function takes output from \code{evian_linear} or \code{evian_logit} as input. It will plot the likelihood intervals for each of the SNPs analyzed. If 1/k interval is significant then it will be colored by the specified color and will remain grey if the interval is not significant.
}


\references{
 Strug, L. J., Hodge, S. E., Chiang, T., Pal, D. K., Corey, P. N., & Rohde, C. (2010). A pure likelihood approach to the analysis of genetic association data: an alternative to Bayesian and frequentist analysis. Eur J Hum Genet, 18(8), 933-941. doi:10.1038/ejhg.2010.47
}
\author{
  Dr. Lisa J Strug \email{lisa.strug@utoronto.ca}
}
\examples{
data(eviandata_linear)
data(evianmap_linear)
\dontshow{
rst1=evian_linear(data=eviandata_linear, bim=evianmap_linear, xcols=10:ncol(eviandata_linear),
ycol=6, covariateCol=c(5,7:9), robust=FALSE, model="additive", m=100, lolim=-0.4,
hilim=0.4, kcutoff = c(32,100), multiThread=1)
}
\donttest{
rst1=evian_linear(data=eviandata_linear, bim=evianmap_linear, xcols=10:ncol(eviandata_linear),
ycol=6, covariateCol=c(5,7:9), robust=FALSE, model="additive",
m=1000, kcutoff = c(32,100), multiThread=1)
}

# Plot the LIs for all 10 SNPs
multiLine_plot(dList=rst1)


}


