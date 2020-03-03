#' The exact distribution of the number of alleles in a m-person DNA mixture
#' 
#' Computes the exact distribution of the number of alleles in a \eqn{m}-person DNA
#' mixture typed with STR loci. For a m-person DNA mixture it is possible to
#' observe \eqn{1,\ldots,2\times m \times L}{1,...,2mL} alleles, where \eqn{L} is 
#' the total number of typed STR loci. The method allows incorporation of the 
#' subpopulation correction, the so-called \eqn{\theta}{theta}-correction, to adjust 
#' for shared ancestry. If needed, the locus-specific probabilities can be obtained using the
#' \code{locuswise} argument.
#' 
#' Computes the exact distribution of the number of alleles for a m-person DNA
#' mixture.
#' 
#' @usage 
#' Pnm_all(m, theta, probs, locuswise = FALSE)
#' Pnm_locus(m, theta, alleleProbs)
#' 
#' @aliases pNoA p.numberofalleles Pnm_locus convolve
#' @param m The number of contributors
#' @param theta The coancestery coefficient
#' @param probs List of vectors with allele probabilities for each locus
#' @param alleleProbs Vectors with allele probabilities
#' @param locuswise Logical. If \code{TRUE} the locus-wise probabilities will be
#' returned. Otherwise, the probability over all loci is returned.
#' @return Returns a vector of probabilities, or a matrix of locuswise
#' probability vectors.
#' @author Torben Tvedebrink, James Curran, Mikkel Andersen
#' @references T. Tvedebrink (2014). 'On the exact distribution of the number of 
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37. 
#' <https://doi.org/10.1007/s00414-013-0951-3>
#' @examples
#' 
#'   ## Simulate some allele frequencies:
#'   freqs <-  structure(replicate(10, { g = rgamma(n = 10, scale = 4, shape = 3); 
#'                                       g/sum(g)
#'                                     },
#'               simplify = FALSE), .Names = paste('locus', 1:10, sep = '.'))
#' 
#'   ## Compute \eqn{\Pr(N(m = 3) = n)}, \eqn{n = 1,\ldots,2 * L *m}, where \eqn{L = 10}
#'   ## here
#'   Pnm_all(m = 2, theta = 0, freqs)
#'   ## Same, but locuswise results
#'   Pnm_all(m = 2, theta = 0, freqs, locuswise = TRUE)
#'   
#' @export
Pnm_all <- function(m, theta, probs, locuswise = FALSE) {
  res = Pnm_all_cpp(m, theta, probs)
  
  if (locuswise == FALSE) {
    res = convolve(res)
  }
  
  return(res)
}

# FIXME
pNoA = Pnm_all
