% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unfazeD.R
\name{gdist}
\alias{gdist}
\title{Genetic Distance}
\usage{
gdist(
  dat,
  afreq,
  coi,
  nmmax,
  nr = 100,
  rval = NULL,
  reval = NULL,
  FUNnm = NULL,
  equalr = FALSE,
  split = "[[:space:][:punct:]]+",
  ...
)
}
\arguments{
\item{dat}{data containing sample genotypes. Various formats allowed.}

\item{afreq}{a list of allele frequencies. Each element of the list
corresponds to a locus.}

\item{coi}{a vector indicating complexity of infection for each sample.}

\item{nmmax}{the maximum number of related strains in a pair of samples.}

\item{nr}{an integer value for the resolution of the grid (\eqn{nr - 1}
values between 0 and 1), over which the likelihood will be calculated.
Ignored if non-null \code{reval} is provided.}

\item{rval}{\eqn{{r}} values for the grid or for evaluation when
\code{equalr} is \code{TRUE}. If \code{NULL}, will be evenly spaced between
0 and 1 and interval \eqn{1/nr}.}

\item{reval}{the grid of \eqn{{r}} combinations, over which the likelihood
will be calculated. A matrix where each column represents a single
combination.}

\item{FUNnm}{potentially a function to select \code{nm} for each pair of
samples.}

\item{equalr}{a logical value. If \code{TRUE}, only equal values of \emph{r} for
different pairs of strains are evaluated.}

\item{split}{the allele separator character string if genotypes are provided
as character strings.}

\item{...}{additional arguments for FUNnm. #*** take out if FUNnm is out}
}
\value{
A lower triangular distance/relatedness matrix in a type of shaped
list. Each pairwise distance contains a scalar, a vector, or a matrix of
\eqn{{r}} estimates. Simple matrix (of type \code{double}) can be returned
if all the elements are of length 1.
}
\description{
Pairwise estimation of relatedness parameters for polyclonal
multiallelic samples.
}
\details{
To be added: data formats that can be processed, formats to be
returned.
}
\seealso{
\code{\link{gdistPair}} for genetic relatedness between two samples
with an option of returning log-likelihood.
}
