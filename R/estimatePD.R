PNdrop <- function(pD, probs, n0, m, loci) {
  i <- 0:((2 * m * loci) - n0)
  ((1 - pD)^n0) * sum(choose(n0 + i, i) * probs[n0 + i] * (pD^i))
}


#' Estimate the drop-out probability based on number of alleles
#' 
#' An inferior may to estimate the drop-out probability compared to using the
#' peak heights from the electropherogram. However, to compare the performance
#' with Gill et al. (2007) this implements a theoretical approach based on
#' their line of arguments.
#' 
#' Computes the \eqn{\Pr(D)}{Pr(D)} that maximises equation (10) in Tvedebrink (2014).
#' 
#' @param n0 Vector of observed allele counts - same length as the number of
#' loci
#' @param m The number of contributors
#' @param pnoa The vector of \eqn{\P(N(m)=n)}{Pr(N(m)=n)} for \eqn{n=1,\ldots,2Lm}, where \eqn{L} is the number
#' of loci and \eqn{m} is the number of contributors OR
#' @param probs List of vectors with allele probabilities for each locus
#' @param theta The coancestery coefficient
#' @param locuswise Logical. Indicating whether computations should be done
#' locuswise.
#' @return Returns the MLE of \eqn{\Pr(D)}{Pr(D)} based on equation (10) in Tvedebrink (2014)
#' @author Torben Tvedebrink
#' @references Gill, P., A. Kirkham, and J. Curran (2007).  LoComatioN: A
#' software tool for the analysis of low copy number DNA profiles.  Forensic
#' Science International 166(2-3): 128 - 138.
#' 
#' T. Tvedebrink (2014). 'On the exact distribution of the number of 
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37. 
#' <https://doi.org/10.1007/s00414-013-0951-3>
#' @examples
#' 
#'   ## Simulate some allele frequencies:
#'   freqs <-  simAlleleFreqs()
#'   ## Assume 15 alleles are observed in a 2-person DNA mixture with 10 loci:
#'   estimatePD(n0 = 15, m = 2, probs = freqs)
#' 
#' @importFrom stats optimise
#' 
#' @export estimatePD
estimatePD <- function(n0, m, pnoa = NULL, probs = NULL, theta = 0, locuswise = FALSE) {
  if (locuswise) {
    pDs <- structure(rep(NA, length(probs)), .Names = names(probs))
    pDs[n0 == 0] <- 1
  } else if (n0 == 0) 
    return(c(`P(D)` = 1))
  if (is.null(pnoa)) {
    if (is.null(probs)) 
      stop("Either 'pNoA' or 'probs' must be provided") else {
      # FIXME:
      pnoa <- Pnm_all(m = m, theta = theta, 
                      probs = probs, locuswise = locuswise)
      #if (locuswise) {
      #  pnoa <- Pnm_locus(m = m, theta = theta, probs)
      #} else {
      #  pnoa <- Pnm_all(m = m, theta = theta, probs)
      #}
    }
  }
  if (locuswise & any(n0 > 2 * m)) 
    stop(paste("One of the number of observed alleles (", (n0[n0 > 2 * m])[1], ") is bigger than the maximum for a ", 
      m, "-person mixture (", 2 * m, ").", sep = "")) else if (!locuswise & any(n0 > length(pnoa))) 
    stop(paste("The number of observed alleles (", n0, ") is bigger than the maximum for a ", 
      m, "-person mixture (", length(pnoa), ").", sep = ""))
  if (locuswise) {
    if (any(n0 != 0)) 
      pDs[n0 != 0] <- mapply(function(x, y) stats::optimise(PNdrop, m = m, probs = x, n0 = y, 
        loci = 1, lower = 0, upper = 1, maximum = TRUE)$maximum, x = pnoa[n0 != 0], 
        y = n0[n0 != 0])
    return(pDs)
  } else {
    loci <- length(pnoa)/(2 * m)
    return(structure(stats::optimise(PNdrop, m = m, probs = pnoa, n0 = n0, loci = loci, lower = 0, 
      upper = 1, maximum = TRUE)$maximum, .Names = c("p(D)")))
  }
}

