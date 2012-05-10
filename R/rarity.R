## Collapses a (m,p)-matrix to a (2*m+p)-vector
mpcollapse <- function(mallele,nloci){
  tmp <- data.frame(m2 = (m <- floor(mallele/2):0),m1=mallele-2*m)
  tmp[tmp$m1<=nloci & tmp$m2<=nloci,]
}

dbCollapse <- function(mpmatrix){
  res <- sapply(1:((nL <- (ncol(mpmatrix)-1))*2+1),function(i) sum(diag(mpmatrix[(MP <- mpcollapse(i-1,nL))$m2+1,,drop=FALSE][,MP$m1+1,drop=FALSE]),na.rm=TRUE))
  names(res) <- 0:(length(res)-1)
  res
}


## Extracts the upper left triangle of a quadratic matrix
up.tri <- function(x,diag=TRUE,droplast=FALSE){
  x <- as.matrix(x)
  res <- (row(x) + col(x)) - 1 <= ncol(x)
  if(!diag) res[(row(x) + col(x)-1) == ncol(x)] <- FALSE
  if(droplast) res[nrow(res),1] <- FALSE
  res
}

## Creates a matrix/list of cell names
dbCats <- function(nloci,vector=FALSE){
  res <- outer(0:nloci,0:nloci,function(i,j) paste(i,j,sep="/"))
  res[!up.tri(res)] <- NA
  if(vector) res <- t(res)[up.tri(res)]
  res
}

## The recursive step of the expected value function. See eq. (2) in Tvedebrink et al.
rare <- function(q){
  S <- ncol(q)
  M <- replicate(S,matrix(0,S+1,S+1),simplify=FALSE)
  M[[1]][1,1] <- q[1,1] # P_{0,1}
  M[[1]][1,2] <- q[2,1] # P_{1,1}
  M[[1]][2,1] <- q[3,1] # P_{2,1}
  for(s in 2:S){
    for(m in 1:(s+1)){
      for(p in 1:(s-m+2)){
        if(m==1 & p>1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[2,s]*M[[s-1]][m,p-1]
        else if(m>1 & p==1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[3,s]*M[[s-1]][m-1,p]
        else if(m==1 & p==1) M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]
        else M[[s]][m,p] <- q[1,s]*M[[s-1]][m,p]+q[2,s]*M[[s-1]][m,p-1]+q[3,s]*M[[s-1]][m-1,p]
      }
    }
  }
  x <- M[[S]]
  dimnames(x) <- list(0:S,0:S)
  x
}

## Computes the P_{0/0}, P_{0/1}, P_{1/0} for a given locus
Ps <- function(p,t,k=rep(0,3)){
  s1 <- sum(p); s2 <- sum(p^2); s3 <- sum(p^3); s4 <- sum(p^4)
  d <- (1+t)*(1+2*t)
  p0 <- t^2*(1-t)*(1-s2) + 2*t*(1-t)^2*(1-2*s2+s3) + (1-t)^3*(1-4*s2+4*s3+2*s2^2-3*s4)
  p1 <- 8*t^2*(1-t)*(1-s2) + 4*t*(1-t)^2*(1-s3) + 4*(1-t)^3*(s2-s3-s2^2+s4)
  p2 <- 6*t^3 + t^2*(1-t)*(2+9*s2) + 2*t*(1-t)^2*(2*s2+s3) + (1-t)^3*(2*s2^2-s4)
  if(all(k==0)) res <- c(p0,p1,p2)/d
  else res <- c(k[3]*p0/d,k[2]*(1-t)*(1-s2)+k[3]*p1/d,k[1]+k[2]*(t+(1-t)*s2)+k[3]*p2/d)
  res
}

## Computes and returns the expected value of the cell counts
dbExpect <- function(probs,theta=0,k=c(0,0,1),n=1,round=FALSE,na=TRUE,vector=FALSE,collapse=FALSE){
  if(length(theta)>1) return(lapply(theta,function(t) dbExpect(probs=probs,theta=t,k=k,n=n,round=round,na=na,vector=vector)))
  if(!is.list(probs) && is.vector(probs)) return(Ps(probs,t=theta,k=k)) ## Previous line handles if more thetas are provided
  ## probs is a list of vectors with each vector being
  ## the allele probabilities for a given locus
  S <- length(probs)
  if(is.character(k)){
    if(!is.na(match(toupper(k),c("UN","AV","FS","FC","PC"))))
      k <- list("UN"=c(0,0,1),"FC"=c(0,1,3)/4,"AV"=c(0,1,1)/2,"PC"=c(0,1,0),"FS"=c(1,2,1)/4)[[toupper(k)]]
    else stop(paste("The value of 'k' (",k,") is not defined - has to be vector c(k2,k1,k0) or string: 'UN' (unrelated), 'FS' (full-siblings), 'AV' (avuncular), 'FC' (first cousins) or 'PC' (parent-child)",sep=""))
  }
  p <- lapply(probs,Ps,t=theta,k=k)
  q <- do.call("cbind",p)
  if(n>1) N <- choose(n,2) else N <- 1
  res <- rare(q)
  if(round) res <- round(res*N)
  else res <- res*N
  if(na) res[!up.tri(res)] <- NA
  if(collapse) return(dbCollapse(res))
  if(vector){
    res <- t(res)[up.tri(res)]
    names(res) <- dbCats(S,vector=TRUE)
  }
  else dimnames(res) <- list(match=0:S,partial=0:S)
  res
}

