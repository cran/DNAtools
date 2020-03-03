#' Compute the posterior probabilities for \eqn{\Pr(m|n_0)}{Pr(m|n0)} for a given prior 
#' \eqn{\Pr(m)}{Pr(m)}.
#' 
#' Compute a matrix of posterior probabilties \eqn{\Pr(m|n_0)}{P(m|n0)} where \eqn{m}
#' ranges from 1 to \eqn{m_{\max}}{\code{m.max}}, and \eqn{n_0}{n0} is
#' \eqn{0,\ldots,2m_{\max}}{0,\ldots,\code{m.max}}. This is done by evaluating
#' \eqn{\Pr(m|n_0)=Pr(n_0|m)Pr(m)/Pr(n)}{Pr(m|n0)=Pr(n0|m)Pr(m)/Pr(n)}, where
#' \eqn{\Pr(n_0|m)}{Pr(n0|m)} is evaluated by \code{\link{pNoA}}.
#' 
#' Computes a matrix of \eqn{\Pr(m|n_0)}{Pr(m|n0)} values for a specific locus.
#' 
#' @param prob Vectors with allele probabilities for the specific locus
#' @param m.prior A vector with prior probabilities (summing to 1), where the
#' length of \code{m.prior} determines the plausible range of \eqn{m}
#' @param m.max Derived from the length of \code{m.prior}, and if
#' \code{m.prior=NULL} a uniform prior is speficied by \code{m.max}:
#' \code{m.prior = rep(1/m.max,m.max)}.
#' @param pnoa.locus A named vector of locus specific probabilities
#' \eqn{P(N(m)=n), n=1,\ldots,2m}{P(N(m)=n), n=1,\ldots,2m}.
#' @param theta The coancestery coefficient
#' @return Returns a matrix \eqn{[\Pr(m|n_0)]}{[Pr(m|n0)]} for 
#' \eqn{m = 1,\ldots,m.max} and \eqn{n_0 = 1,\ldots,2m.max}{n_0 = 1,\ldots,2m.max}.
#' @author Torben Tvedebrink, James Curran
#' @references T. Tvedebrink (2014). 'On the exact distribution of the number of 
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37. 
#' <https://doi.org/10.1007/s00414-013-0951-3>
#' @examples
#' 
#'   ## Simulate some allele frequencies:
#'   freqs <-  simAlleleFreqs()
#'   
#'   ## Compute Pr(m|n0) for m = 1, ..., 5 and n0 = 1, ..., 10 for the first locus:
#'   pContrib_locus(prob = freqs[[1]], m.max = 5)
#' 
#' @export pContrib_locus
pContrib_locus <- function(prob = NULL, m.prior = NULL, m.max = 8, 
                           pnoa.locus = NULL, theta = 0) {
  if (is.null(m.prior)) 
    m.prior <- rep(1/m.max, m.max)
  
  if (length(m.prior) != m.max) 
    m.max <- length(m.prior)
  
  if (!is.null(pnoa.locus)) {
    if (!is.list(pnoa.locus)) 
      stop("'pnoa.locus' must be a named list of locus probabilities, e.g. list('1'=locus.m1,'2'=locus.m2)")
    
    if (!all(1:m.max %in% names(pnoa.locus)) & is.null(prob)) 
      stop("If not all locuswise probabilities are supplied for 1,...,m.max, then a vector of allele frequencies must be supplied.")
    
    PNOAlocus <- structure(replicate(m.max, NULL, simplify = FALSE), .Names = 1:m.max)
    for (i in 1:m.max) {
      if (i %in% names(pnoa.locus)){ 
        PNOAlocus[[i]] <- c(pnoa.locus[[paste(i)]], rep(0, 2 * (m.max - i)))
      }else{
        PNOAlocus[[i]] <- c(Pnm_locus(m = i, theta = theta, alleleProbs = prob), 
                            rep(0, 2 * (m.max - i)))
      }
    }
  } else {
    PNOAlocus <- lapply(1:m.max, function(i) c(Pnm_locus(m = i, theta = theta, 
                                                         alleleProbs = prob), 
                                                         rep(0, 2 * (m.max - i))))
  }
  PNOAlocus <- do.call("rbind", PNOAlocus)
  ## Compute the posterior: P(m|n) = P(n|m)*P(m)/P(n); where the computations are done for
  ## n=1,...,2*m.max The denominator in vector format: P(n) = (P(n=1),..,P(n=2*m.max)) =
  ## (sum{j=1}^m.max P(n=1|m=j)*P(m=j), ..., sum{j=1}^m.max P(n=2*m.max|m=j)*P(m=j))
  denominator <- as.numeric(m.prior %*% PNOAlocus)
  ## Compute numerator: a vector of vectors (matrix): [P(n=1|m=1)*p(m=1) ... P(n=1|m=M)*P(m=M);
  ## .... ; P(n=2M|m=1)*p(m=1) ... P(n=2M|m=M)*P(m=M)]
  numerator <- t(PNOAlocus) * m.prior
  ## Return a matrix of posterior probabilities
  structure(numerator/denominator, .Dimnames = list(n0 = 1:(2 * m.max), `P(m|n0)\n   m` = 1:m.max))
}

### Computes the posterior probabilities of m given a specific vector of observed locus counts


#' Compute the posterior probabilities for P(m|n0) for a given prior P(m) and
#' observed vector n0 of locus counts
#' 
#' where m ranges from 1 to \eqn{m_{\max}}{\code{m.max}} and \eqn{n_0}{n0} is
#' the observed locus counts.
#' 
#' Computes a vector P(m|n0) evaluated over the plausible range 1,...,m.max.
#' 
#' @param n0 Vector of observed allele counts - same length as the number of
#' loci.
#' @param probs List of vectors with allele probabilities for each locus
#' @param m.prior A vector with prior probabilities (summing to 1), where the
#' length of \code{m.prior} determines the plausible range of m
#' @param m.max Derived from the length of \code{m.prior}, and if
#' \code{m.prior=NULL} a uniform prior is speficied by \code{m.max}:
#' \code{m.prior = rep(1/m.max,m.max)}.
#' @param theta The coancestery coefficient
#' @return Returns a vector P(m|n0) for m=1,...,m.max
#' @author Torben Tvedebrink, James Curran
#' @references T. Tvedebrink (2014). 'On the exact distribution of the number of 
#' alleles in DNA mixtures', International Journal of Legal Medicine; 128(3):427--37. 
#' <https://doi.org/10.1007/s00414-013-0951-3>
#' @examples
#' 
#'   ## Simulate some allele frequencies:
#'   freqs <-  simAlleleFreqs()
#'   m <- 2
#'   n0 <- sapply(freqs, function(px){
#'                             peaks = unique(sample(length(px),
#'                                              size = 2 * m,
#'                                              replace = TRUE,
#'                                              prob = px))
#'                             return(length(peaks))
#'                        })
#'   ## Compute P(m|n0) for m=1,...,4 and the sampled n0
#'   pContrib(n0=n0,probs=freqs,m.max=4)
#' 
#' @export pContrib
pContrib <- function(n0, probs = NULL, 
                     m.prior = rep(1/m.max, m.max), m.max = 8, 
                     theta = 0) {
  
  ## Check that the length of the prior is the same as m.max,
  ## if it isn't then change m.max to the length of the prior.
  ## JMC: deleted the if statement, because if it matches then there will
  ## be no difference
  m.max = length(m.prior)
  
  ## Calculate the probability tables for Pr(n0 | m = j)
  ## n0 is a vector with the same length as the number of loci
  ## therefore Pr(n0 | m = j) = \prod_{l = 1}^{L}Pr(n_{0,l} | m = j)
  
  pr_n0_m = lapply(1:m.max, function(j){
    Pnm_all(m = j, theta = theta, probs = probs, locuswise = TRUE)
  })
  
  ## The previous step gives a list with m.max elements labelled m = 1,...,m.max
  ## with the probability of seeing 1,.., 2 * m peaks at each locus
  ## so now we need to extract the exact value of Pr(n0 | m = j), which isn't hard
  ## but we need to make sure it is zero when n0 > 2 * m
  ## I am going to use a loop here, because m.max is typically small
  
  nLoci = length(probs)
  
  jointProb = rep(0, m.max)
  
  for(j in 1:m.max){
    locusProb = rep(0, nLoci)
    
    for(l in 1:nLoci){
      if(n0[l] <= 2 * j){
        locusProb[l] = pr_n0_m[[j]][l, n0[l]]
      }
    }
    
    jointProb[j] = prod(locusProb) * m.prior[j]
  }
  
  
  numerator <- jointProb ## holds P(m|n0)*P(m) for m=1,...,length(m.prior)
  ## The denominator is simply the sum over m in the numerator:
  structure(numerator/sum(numerator), .Names = paste(paste("P(m=", 1:m.max, sep = ""), "|n0)", 
    sep = ""))
}

