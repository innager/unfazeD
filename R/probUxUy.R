#' Calculate likelihood for a pair \ifelse{html}{\out{U<sub>x</sub>}}{\eqn{U_x}}, \ifelse{html}{\out{U<sub>y</sub>}}{\eqn{U_y}}
#'
#' Calculates log-likelihood for a pair of samples at a single locus.
#'
#' @details \code{probUxUyEqr} is used for a special case when all pairs of
#'   strains are assumed to be related with the same parameter \eqn{r}. For that
#'   case, fast calculations are performed, and the results are equivalent to
#'   the general case with each value of \code{rval} is repeated \code{nm}
#'   times.
#'
#' @param Ux,Uy sets of unique alleles for two samples at a given locus. Vectors
#'   of indices corresponding to ordered probabilities in \code{probs}.
#' @param nx,ny complexity of infection for two samples. Vectors of length 1.
#' @param probs a vector of population allele frequencies at a given locus.
#'   (Should sum to one - not checked!)
#' @param reval a matrix containing the grid of \eqn{{r}} combinations, over
#'   which the likelihood is to be calculated. Each column represents a single
#'   combination.
#'
#' @return A vector of log-likelihood values for each evaluated \eqn{{r}}
#'   combination (or a value for \code{probUxUyEqr}.
#'
#' @export

probUxUy <- function(Ux, Uy, nx, ny, probs, reval) {
  nm <- nrow(reval)
  if (nm > (min(nx, ny))) {
    stop("number of related strains greater than min(nx, ny)")
  }
  ixy <- which(Ux %in% Uy)
  iyx <- which(Uy %in% Ux)
  logj   <- log(1:max(nx, ny))         # starts with 1
  factj  <- lgamma(0:max(nx, ny) + 1)  # starts with 0

  return(.Call("llik", as.integer(Ux), as.integer(Uy), as.integer(ixy),
               as.integer(iyx), as.integer(nx), as.integer(ny),
               as.double(probs), as.double(logj), as.double(factj),
               as.double(reval), as.integer(nm), PACKAGE = "unfazeD"))
}

#' @rdname probUxUy
#' @param nm   the number of related strains.
#' @param rval a vector of \eqn{r} values, over which the likelihood is to be
#'   calculated.
#' @export
#'
probUxUyEqr <- function(Ux, Uy, nx, ny, probs, rval, nm) {
  if (nm > (min(nx, ny))) {
    stop("number of related strains greater than min(nx, ny)")
  }
  ixy <- which(Ux %in% Uy)
  iyx <- which(Uy %in% Ux)
  logj   <- log(1:max(nx, ny))         # starts with 1
  factj  <- lgamma(0:max(nx, ny) + 1)  # starts with 0

  return(.Call("llikEqr", as.integer(Ux), as.integer(Uy),
               as.integer(ixy), as.integer(iyx), as.integer(nx), as.integer(ny),
               as.double(probs), as.double(logj), as.double(factj),
               as.double(rval), as.integer(nm), PACKAGE = "unfazeD"))
}
