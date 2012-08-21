old.dbCompare <- function(x,profile=NULL,hit=7,trace=TRUE,vector=FALSE,collapse=FALSE,threads=2,compareF=FALSE,imposeF=FALSE){
  for(i in 1:ncol(x)) x[,i] <- paste(x[,i])
  nl <- (ncol(x)-1)/2
  if(imposeF){
    compareF <- TRUE
    for(j in (1:nl)*2+1) x[x[,j-1]==x[,j],j] <- "0" ## Force all homozygotes to have aF genotype (here '0' place the role of F)
  } 
  if(!is.null(profile)){
    nx <- names(x)
    ## The length of the profile is odd, no id given:
    if((length(profile)%%2)==0) profile <- c("profile",paste(profile))
    names(profile) <- nx
    x <- rbind(paste(profile),x)
    single <- 1
  }
  else single <- 0
  if(trace) progress <- 1 else progress <- 0
  ids <- unlist(x[,1])
  nid <- names(x)[1]
  x <- do.call("paste",c(x=x,sep="\t"))
  if(compareF){
    if(threads==1) res <- .Call("fcompare", x, c(as.integer(nl),as.integer(hit),as.integer(progress),as.integer(single)), PACKAGE = "DNAtools")
    else res <- .Call("fmcompare", x, c(as.integer(nl),as.integer(hit),as.integer(progress),as.integer(single),as.integer(threads)), PACKAGE = "DNAtools")
    allCats <- as.vector(t(outer(0:nl,0:nl,FUN=paste,sep="/")))
    remCats <- as.vector(t(outer(0:nl,0:nl,FUN=function(x,y,n=nl) (x+y)<=n)))
    result <- list(m=matrix(res$M,(nl+1)*(nl+1),(nl+1)*(nl+1),byrow=FALSE,dimnames=list(Genuine=allCats,Wildcard=allCats))[remCats,remCats])
  }
  else{
    if(threads==1) res <- .Call("compare", x, c(as.integer(nl),as.integer(hit),as.integer(progress),as.integer(single)), PACKAGE = "DNAtools")
    else res <- .Call("mcompare", x, c(as.integer(nl),as.integer(hit),as.integer(progress),as.integer(single),as.integer(threads)), PACKAGE = "DNAtools")
    result <- list(m=matrix(res$M,nl+1,nl+1,byrow=TRUE,dimnames=list(match=0:nl,partial=0:nl)))
  }
  call <- list(loci=nl,single=single,collapse=collapse)
  if(length(res$row1)>0){
    if(compareF) result <- c(result,list(hits=data.frame(id1=ids[res$row1],id2=ids[res$row2],match=res$matches,partial=res$partial,Fmatch=res$fmatches,Fpartial=res$fpartial)))
    else result <- c(result,list(hits=data.frame(id1=ids[res$row1],id2=ids[res$row2],match=res$matches,partial=res$partial)))
    names(result$hits)[1:2] <- paste(nid,1:2,sep="")
    call <- c(call,list(hit=hit))
  }
  if(collapse){
    if(compareF){
      result$m <- varCollapse(result$m,nl=nl)
      dimnames(result$m) <- list(Genuine=1:nrow(result$m)-1,Wildcard=1:ncol(result$m)-1)
      if(vector) result$m <- wildCollapse(result$m)
    }
    else{
      result$m <- dbCollapse(result$m)
      names(result$m) <- 0:(length(result$m)-1)
    }
  }
  if(vector & !collapse){
    if(compareF) result$m <- result$m[up.tri(result$m)]
    else result$m <- t(result$m)[up.tri(result$m)]
  }
  attributes(result)$call <- call
  attributes(result)$class <- "dbcompare"
  result
}

print.old.dbcompare <- function(x,...){
  if(is.matrix(x$m)) x$m[!up.tri(x$m)] <- NA
  if(length(attributes(x)$names)==1){
    if(is.matrix(x$m)) print.table(x$m,...)
    else print.default(x$m,...)
  }
  else{
    if(is.matrix(x$m)){
      cat("Summary matrix\n")
      print.table(x$m,...)
    }
    else{
      cat("Summary vector\n")
      print.default(x$m,...)
    }
    cat(paste("\nProfiles with at least",attributes(x)$call$hit,"matching loci\n"))
    x$hits <- x$hits[order(x$hits$match,x$hits$partial,decreasing=TRUE),]
    rownames(x$hits) <- 1:nrow(x$hits)
    print.data.frame(x$hits,...)
  }
}

