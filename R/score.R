score = function(prof1, prof2, bWildCard = FALSE, bRareAllele = FALSE, nRareCode = 99) {
  ## assumes that profiles have id as the first entry and then 2*n alleles
  nLoci = (length(prof1) - 1)/2
  
  i1 = which(!(nRareCode %in% prof1[(2 * (1:nLoci)):(2 * (1:nLoci) + 1)]) && (prof1[2 * (1:nLoci)] > 
    prof1[2 * (1:nLoci) + 1]))
  if (any(i1)) {
    tmp = prof1[2 * i1]
    prof1[2 * i1] = prof1[2 * i1 + 1]
    prof1[2 * i1 + 1] = tmp
  }
  
  i2 = which(!(nRareCode %in% prof2[(2 * (1:nLoci)):(2 * (1:nLoci) + 1)]) && (prof2[2 * (1:nLoci)] > 
    prof2[2 * (1:nLoci) + 1]))
  if (any(i2)) {
    tmp = prof2[2 * i2]
    prof2[2 * i2] = prof2[2 * i2 + 1]
    prof2[2 * i2 + 1] = tmp
  }
  
  # mikl scores = .Call('score_cpp')
  scores = score_rcpp(prof1 = as.integer(prof1[-1] * 10), prof2 = as.integer(prof2[-1] * 10), 
    numLoci = nLoci, useWildCard = bWildCard, useRareAllele = bRareAllele)
  
  g1 = paste(prof1[2 * (1:nLoci)], prof1[2 * (1:nLoci) + 1], sep = "/")
  g2 = paste(prof2[2 * (1:nLoci)], prof2[2 * (1:nLoci) + 1], sep = "/")
  
  g1 = gsub(as.character(nRareCode), "R", g1)
  g2 = gsub(as.character(nRareCode), "R", g2)
  
  return(data.frame(g1 = g1, g2 = g2, score = scores))
}
