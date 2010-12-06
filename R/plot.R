plot.dbcompare <- function(x,log="y",las=3,xlab="Match/Partial",ylab="Counts",...){
  nl <- attributes(x)$call$loci
  levs <- dbCats(nl,vector=TRUE)
  if(is.matrix(x$m)) mvec <- t(x$m)[up.tri(x$m)]
  else mvec <- x$m
  mvec[mvec==0] <- NA  
  plot(1:length(levs),mvec,axes=FALSE,xlab=xlab,ylab=ylab,log=log,...)
  axis(1,at=1:length(levs),labels=levs,las=las)
  axis(2)
  box()
}
