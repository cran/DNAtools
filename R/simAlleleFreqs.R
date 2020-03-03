#' Simulate Allele Frequencies
#' 
#' Simulate some allele frequencies using Dirichlet Random variables
#'
#' @param nLoci \eqn{L} the number of loci in the multiplex
#' @param allelesPerLocus  the number of alleles per locus
#' @param shape the shape parameter
#'
#' @return a list with elements \code{locus.}\eqn{l} where \eqn{l=1,\ldots,L}, each 
#' of which are vectors of length \code{allelesPerLocus[l]}, consisting of allele 
#' frequencies for that locus
#' @export
#'
#' @examples
#' set.seed(123)
#' simAlleleFreqs()
#' 
#' @importFrom stats rgamma
simAlleleFreqs = function(nLoci = 10, allelesPerLocus = rep(10, nLoci),
                          shape = rep(3, nLoci)){
  
  if(length(allelesPerLocus) != nLoci){
    stop(paste0("If the number of alleles at each locus is different,\n",
                "then the number of alleles should be specified for each locus."))
  }
  
  if(length(shape) != nLoci){
    stop(paste0("If the shape parameters are different for each locus,\n",
                "then a shape parameter should be specified for each locus."))
  }
  
  freqs = lapply(allelesPerLocus, function(nA){
                 stats::rgamma(nA, shape = shape, rate = 1)
  })
    
  freqs = lapply(freqs, function(l){l/sum(l)})
  names(freqs) = paste("locus", 1:length(freqs), sep = ".")

  return(freqs)
}
