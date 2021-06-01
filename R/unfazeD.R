#' Genetic Distance for a Pair of Samples
#'
#' Calculates likelihood over a range of values for relatedness parameter
#' \eqn{{r}} and provides an estimate for a pair of samples.
#'
#' @param pair   a list of length two containing two samples.
#' @param afreq  a list of allele frequencies. Each element of the list
#'   corresponds to a locus.
#' @param coi    a vector indicating complexity of infection for each sample.
#' @param nr     an integer value for the resolution of the grid (\eqn{nr - 1}
#'   values between 0 and 1), over which the likelihood will be calculated.
#'   Ignored if non-null \code{reval} is provided.
#' @param nm     the number of related pairs of strains for evaluation.
#' @param rval  \eqn{{r}} values for the grid or for evaluation when
#'   \code{equalr} is \code{TRUE}. If \code{NULL}, will be evenly spaced between
#'   0 and 1 and interval \eqn{1/nr}.
#' @param reval  the grid of \eqn{{r}} combinations, over which the likelihood
#'   will be calculated. A matrix where each column represents a single
#'   combination.
#' @param equalr a logical value. If \code{TRUE}, only equal values of _r_ for
#'   different pairs of strains are evaluated.
#' @param out    a character string for the type of results to be returned. If
#'   \code{"mle"}, an estimate is returned, if \code{"llik"} - a vector of
#'   log-likelihood for each combination.
#'
#' @return If \code{out} is \code{"mle"}, a vector of length 1 if \code{equalr}
#'   is \code{TRUE} or length \code{nm} otherwise, containing estimated
#'   \eqn{{r}} (or vector/matrix if not just first). If \code{out} is
#'   \code{"llike"}, a vector of log-likelihood values for each evaluated
#'   combination.
#'
#' @seealso \code{\link{gdist}} for processing multi-sample data in various
#'   formats and estimating pairwise relatedness.
#' @export
#' @useDynLib unfazeD

#*** add NA handling
gdistPair <- function(pair, afreq, coi, nr = 1e2, nm = min(coi), rval = NULL,
                      reval = NULL, equalr = FALSE, out = "mle") {  # "mle", "llik", "all"
  nloc <- length(afreq)

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateReval(nm, rval)[[nm]]
    }
    neval <- ncol(reval)
  }

  llik <- rep(0, neval)
  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))
    if (length(Ux) == 0 ||                   # NA [or all 0's - not in true]
        length(Uy) == 0) {#||                # NA [or all 0's - not in true]
#        any(afreq[[t]][unique(c(Ux, Uy))] <= 0)) {  # not in true genotypes
      next                                   # no data (NA)
    }
    if (equalr) {
      llikt <- probUxUyEqr(Ux, Uy, coi[1], coi[2], afreq[[t]], rval, nm)
    } else {
      llikt <- probUxUy(   Ux, Uy, coi[1], coi[2], afreq[[t]], reval)
    }
    if (llikt[1] == -Inf) {
      if (all(llikt == -Inf)) {
        next
      } else {
        iinf <- which(llikt == -Inf)
        writeLines(paste("\n0 prob for r combinations",  # not added if 1st!!
                         paste(iinf, collapse = " ")))
      }
    } else {
      llik <- llik + llikt
    }
  }

  if (out == "llik") {
    return(llik)
  } else {
    imax <- which(llik == max(llik))
    if (equalr) {
      return(mean(rval[imax]))  # [1] first only; reval[imax] - all; mean()
    } else {
      return(rowMeans(reval[, imax, drop = FALSE]))
      # reval[, imax[1]]
    }
  }
}

#' Genetic distance for data with genotyping errors
#'
#' Distance when estimated COI is less than |Ux|: COI increased to |Ux| for a
#' given locus only. ***Handling of the situation where alleles with estimated
#' allele frequency of 0 are present: probUxUy crashes only if present in both
#' Ux and Uy; however, if present in either one, the likelihood is 0, so no need
#' to calculate.
#'
#' @param alpha significance level.
#' @return
#'   * If \code{out = "mle"}: a vector of length 1 if \code{equalr = TRUE}
#'   or of length \code{nm} otherwise, containing estimated \eqn{{r}} (or vector
#'   /matrix if not just first).
#'   * If \code{out = "llike"}, a vector of log-likelihood values for each
#'   evaluated combination.
#'   * If \code{out = "all"}, a list containing an estimate, log-likelihood,
#'   maximum log-likelihood, \eqn{{r}} values corresponding to the acceptance
#'   region determined by the significance level, and the size of that region
#'   (as a proportion of all evaluated).
#' @inherit gdistPair params
#' @export

# out options: "mle", "llik", "all" (*** gdistPair1 only - add to others?)
gdistPair1 <- function(pair, afreq, coi, nr = 1e2, nm = min(coi), rval = NULL,
                       reval = NULL, equalr = FALSE, out = "mle",
                       alpha = 0.05) {
  nloc <- length(afreq)

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateReval(nm, rval)[[nm]]
    }
    neval <- ncol(reval)
  }

  llik <- rep(0, neval)
  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))
    #*** shouldn't happen with our error simulation
    if (length(Ux) == 0 ||                   # NA or all 0's
        length(Uy) == 0 ||                   # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] == 0)) {  # likelihood = 0
      next                                   # likelihood = 0 or no data
    }
    coix <- max(coi[1], length(Ux))
    coiy <- max(coi[2], length(Uy))
    if (equalr) {
      llikt <- probUxUyEqr(Ux, Uy, coix, coiy, afreq[[t]], rval, nm)
    } else {
      llikt <- probUxUy(   Ux, Uy, coix, coiy, afreq[[t]], reval)
    }
    #*** check and maybe remove
    if (llikt[1] == -Inf) {
      if (all(llikt == -Inf)) {
        next
      } else {
        iinf <- which(llikt == -Inf)
        writeLines(paste("\n0 prob for r combinations",  # not added if 1st!!
                         paste(iinf, collapse = " ")))
      }
    } else {
      llik <- llik + llikt
    }
  }

  if (tolower(out) == "llik") {
    return(llik)
  }
  imax <- which(llik == max(llik))
  if (equalr) {
    est <- mean(rval[imax])  # [1] first only; reval[imax] - all; mean()
  } else {
    est <- rowMeans(reval[, imax, drop = FALSE])
    # reval[, imax[1]]
  }
  if (tolower(out == "mle")) {
    return(est)
  }

  qchi <- stats::qchisq(1 - alpha, df = ifelse(equalr, 1, nrow(reval)))
  cutoff  <- max(llik) - qchi/2
  itop    <- llik >= cutoff
  proptop <- sum(itop)/neval
  if (equalr) {
    rtop <- rval[itop]
  } else {
    rtop <- reval[, itop, drop = FALSE]
  }
  return(list(mle = est, llik = llik, maxllik = llik[imax[1]], rtop = rtop,
              proptop = proptop))
}

#' Genetic distance for data with genotyping errors
#'
#' Distance when estimated COI is less than |Ux|: likelihoods for all subsets of
#' size = COI of Ux are calculated; geometric mean taken.
#'
#' @inherit gdistPair return params
#' @export
#'

gdistPair2 <- function(pair, afreq, coi, nr = 1e2, nm = min(coi), rval = NULL,
                       reval = NULL, equalr = FALSE, out = "mle") {  # "mle", "llik"
  nloc <- length(afreq)

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateReval(nm, rval)[[nm]]
    }
    neval <- ncol(reval)
  }

  llik <- rep(0, neval)
  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))  # sorted
    if (length(Ux) == 0 ||                   # NA or all 0's
        length(Uy) == 0 ||                   # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] <= 0)) {  #*** < check 0 separately?
      next                                   # likelihood = 0 or no data
    }
    Uxcomb <- getComb(Ux, coi[1])
    Uycomb <- getComb(Uy, coi[2])
    ncx <- ncol(Uxcomb)
    ncy <- ncol(Uycomb)
    llikt <- matrix(0, ncx*ncy, neval)
    icomb <- 1
    for (icx in 1:ncx) {
      for (icy in 1:ncy) {
        if (equalr) {
          llikt[icomb, ] <- probUxUyEqr(Uxcomb[, icx], Uycomb[, icy],
                                         coi[1], coi[2], afreq[[t]], rval, nm)
        } else {
          llikt[icomb, ] <- probUxUy(   Uxcomb[, icx], Uycomb[, icy],
                                        coi[1], coi[2], afreq[[t]], reval)
        }
        icomb <- icomb + 1
      }
    }
    llikt <- colMeans(llikt, na.rm = TRUE)  # or FALSE?
    if (llikt[1] == -Inf) {
      if (all(llikt == -Inf)) {
        next
      } else {
        iinf <- which(llikt == -Inf)
        writeLines(paste("\n0 prob for r combinations",  # not added if 1st!!
                         paste(iinf, collapse = " ")))
      }
    } else {
      llik <- llik + llikt
    }
  }

  if (out == "llik") {
    return(llik)
  } else {
    imax <- which(llik == max(llik))
    if (equalr) {
      return(mean(rval[imax]))  # [1] first only; reval[imax] - all; mean()
    } else {
      return(rowMeans(reval[, imax, drop = FALSE]))
      # reval[, imax[1]]
    }
  }
}

#' Genetic Distance
#'
#' Pairwise estimation of relatedness parameters for polyclonal multiallelic
#' samples.
#'
#' @details To be added: data formats that can be processed, formats to be
#'   returned.
#'
#' @param dat data containing sample genotypes. Various formats allowed.
#' @param nmmax the maximum number of related strains in a pair of samples.
#' @param split the allele separator character string if genotypes are provided
#'   as character strings.
#' @param FUNnm potentially a function to select \code{nm} for each pair of
#'   samples.
#' @param ...   additional arguments for FUNnm. #*** take out if FUNnm is out
#' @inheritParams gdistPair
#'
#' @return A lower triangular distance/relatedness matrix in a type of shaped
#'   list. Each pairwise distance contains a scalar, a vector, or a matrix of
#'   \eqn{{r}} estimates. Simple matrix (of type \code{double}) can be returned
#'   if all the elements are of length 1.
#'
#' @seealso \code{\link{gdistPair}} for genetic relatedness between two samples
#'   with an option of returning log-likelihood.
#' @export

#*** maybe add an option where we want a function for nm, not just a nmmax
#    e.g. max(coi) - 1; provide a function?
gdist <- function(dat, afreq, coi, nmmax, nr = 1e2, rval = NULL, reval = NULL,
                  FUNnm = NULL,
                  equalr = FALSE, split = "[[:space:][:punct:]]+", ...) {
  if (inherits(dat, "matrix")) {
    nsmp <- nrow(dat)
    nloc <- ncol(dat)
    if (typeof(dat) == "list") {
      format <- "matlist"
    } else {
      format <- "matstr"
    }
  } else {
    nsmp <- length(dat)
    nloc <- length(dat[[1]])
    format <- "listlist"
  }

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateReval(1:nmmax, rval)
    }
  }
  res  <- matrix(list(NA), nsmp, nsmp)

  for (ix in 1:nsmp) {
    for (iy in 1:(ix)) {
      nm <- min(coi[ix], coi[iy], nmmax)       # or FUNnm() - or AFTER ix == iy
      if (ix == iy) {                          # diagonal
        res[[ix, ix]] <- rep(1, ifelse(equalr, 1, min(coi[ix], nmmax)))
        next
      }
      llik <- rep(0, ifelse(equalr, neval, ncol(reval[[nm]])))
      for (t in 1:nloc) {
        if (format == "matstr") {
          alleles <- names(afreq[[t]])
          if (is.null(alleles))
            stop("Please provide names for allele frequencies")
          Ux  <- unlist(strsplit(dat[ix, t], split = split))
          Ux  <- sort(match(Ux, alleles))
          Uy  <- unlist(strsplit(dat[iy, t], split = split))
          Uy  <- sort(match(Uy, alleles))
        } else if (format == "matlist") {
          Ux <- which(as.logical(dat[[ix, t]]))
          Uy <- which(as.logical(dat[[iy, t]]))
        } else if (format == "listlist") {
          Ux <- which(as.logical(dat[[ix]][[t]]))
          Uy <- which(as.logical(dat[[iy]][[t]]))
        }
        if (length(Ux) == 0 || is.na(Ux[1]) ||   # NA or all 0's
            length(Uy) == 0 || is.na(Uy[1]) ||   # NA or all 0's
            any(afreq[[t]][unique(c(Ux, Uy))] <= 0)) {  #*** < check 0 separately?
          next                                   # likelihood = 0 or no data
        }
        if (equalr) {
          llikt <- probUxUyEqr(Ux, Uy, coi[ix], coi[iy], afreq[[t]], rval, nm)
        } else {
          llikt <- probUxUy(Ux, Uy, coi[ix], coi[iy], afreq[[t]], reval[[nm]])
        }
        if (any(llikt == -Inf)) {
          writeLines(paste("samples", ix, "and", iy, "- 0 prob for some r
                           combinations at locus", t))
          if (all(llikt == -Inf)) {       # -Inf is added otherwise
            next
          }
          llik <- llik + llikt            # lik = 0 if -Inf at any locus
        }
      }
      imax <- which(llik == max(llik))
      if (equalr) {
        res[[ix, iy]] <- mean(rval[imax])  # [1] first; rval[imax] all; mean()
      } else {
        res[[ix, iy]] <- rowMeans(reval[[nm]][, imax, drop = FALSE])
        # reval[[nm]][, imax[1]] first
      }
    }
  }
  if (all(sapply(res, length) == 1)) {
    res <- matrix(unlist(res), nrow(res), ncol(res))  # return simple matrix
  }
  return(res)
}

#' Generate a grid of parameter values to evaluate over
#'
#' @param Mmax an integer, maximum \code{M}.
#' @param rval  \eqn{{r}} values for the grid.
#' @param nr   an integer vector of length \code{Mmax} indicating the fineness
#'   of the grid.
#' @return A list of length \code{max(Meval)}. Indices of the elements
#'   correspond to values \code{i} in \code{Meval}; each such element is a
#'   matrix with \code{i} rows. Other elements (if any) are \code{NULL} except
#'   for the first, which is always included.
#'
#' @export
# list returned, elements correspond to Meval; M = 1 included always
# if rval and nr are both provided, rval used for M = 1, nr used for M > 1
generateReval <- function(Mmax, rval = NULL, nr = NULL) {
  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr[1]), ceiling(log(nr[1], 10)))
  }
  reval <- list(matrix(rval, 1))
  for (M in 2:Mmax) {
    if (!is.null(nr)) {
      rval <- round(seq(0, 1, 1/nr[M]), ceiling(log(nr[M], 10)))
    }
    revalm <- as.matrix(expand.grid(rep(list(rval), M)))
    for (k in 1:nrow(revalm)) {         # faster than apply()
      revalm[k, ] <- sort(revalm[k, ])
    }
    reval[[M]] <- t(unique(revalm))
  }
  return(reval)
}

# get combinations of alleles when length(Ux) > coix
getComb <- function(Ux, coix) {
  if (length(Ux) <= coix) {
    return(matrix(Ux, ncol = 1))
  }
  return(utils::combn(Ux, coix))
}


