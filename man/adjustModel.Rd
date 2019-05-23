\name{adjustModel}
\alias{adjustModel}
\title{Genotype coding adjusement}
\description{
 This is a helper function that adjusts the genotye coding scheme based on the genetic model specified.
}
\usage{
  adjustModel(data_nomiss,model)}
\arguments{
  \item{data_nomiss}{a data frame that contains the phenotype, covariates and genotype columns. Genotypes need to be coded as 0/1/2.}
  \item{model}{The model specified. Must be a string from one of the following: \code{additive}, \code{dominant}, \code{recessive}, or \code{overdominance}}
}
\details{
  \code{adjustModel} is an interior function that adjusts the genotype coding based on the genetic model specified. The coding scheme for different genetic models can be found in \code{\link{calculateEvianMLE}}.
}
\value{
This function returns the data frame with the same columns but changed genotype coding based on the genetic model specified.
}

\keyword{models}

