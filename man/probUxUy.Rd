% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probUxUy.R
\name{probUxUy}
\alias{probUxUy}
\alias{probUxUyEqr}
\title{Calculate likelihood for a pair \ifelse{html}{\out{U<sub>x</sub>}}{\eqn{U_x}}, \ifelse{html}{\out{U<sub>y</sub>}}{\eqn{U_y}}}
\usage{
probUxUy(Ux, Uy, nx, ny, probs, reval)

probUxUyEqr(Ux, Uy, nx, ny, probs, rval, nm)
}
\arguments{
\item{Ux, Uy}{sets of unique alleles for two samples at a given locus. Vectors
of indices corresponding to ordered probabilities in \code{probs}.}

\item{nx, ny}{complexity of infection for two samples. Vectors of length 1.}

\item{probs}{a vector of population allele frequencies at a given locus.
(Should sum to one - not checked!)}

\item{reval}{a matrix containing the grid of \eqn{{r}} combinations, over
which the likelihood is to be calculated. Each column represents a single
combination.}

\item{rval}{a vector of \eqn{r} values, over which the likelihood is to be
calculated.}

\item{nm}{the number of related strains.}
}
\value{
A vector of log-likelihood values for each evaluated \eqn{{r}}
combination (or a value for \code{probUxUyEqr}.
}
\description{
Calculates log-likelihood for a pair of samples at a single locus.
}
\details{
\code{probUxUyEqr} is used for a special case when all pairs of
strains are assumed to be related with the same parameter \eqn{r}. For that
case, fast calculations are performed, and the results are equivalent to
the general case with each value of \code{rval} is repeated \code{nm}
times.
}
