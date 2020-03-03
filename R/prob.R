prob = function(combinations, freqs, rare, threshold, theta = 0) {
  bIsAll = any(grepl("all", combinations, ignore.case = TRUE))
  
  if (bIsAll & length(combinations) > 1) {
    warning("All option used, all other combinations will be ignored")
    combinations = "all"
  }
  
  res <- Prob(vstrCombs = combinations, q = freqs, R = rare, r = threshold, t = theta)
  return(res)
  # .Call('_DNAtools2_Prob', combinations, freqs, rare, threshold, theta, PACKAGE =
  # 'DNAtools2')
}
