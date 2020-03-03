#' Compare DNA profiles
#' 
#' Compare DNA profiles
#' 
#' Computes the distance between DNA profiles in terms of matching and
#' partially-matching STR loci.
#' 
#' @param x Database with DNA profiles. The database format is expected to be a
#' data frame with each column containing an allelic number such that for each
#' DNA marker there are two columns in the data frame. See
#' \code{data(dbExample)} for an example of the format.
#' @param profiles One or more profiles to be compared with all profiles in the
#' database. Input is a vector, matrix or data frame of same length/width as a
#' row in the database \code{x}.  If profiles is non-null only one CPU will be
#' used. In case threads>1 a warning will be given but computations performed
#' using single core.
#' @param hit The number of matching loci for further investigation
#' @param trace Shows a progress bar
#' @param vector Logical. Whether the result should be returned as vector or a
#' matrix. Note if 'collapse' is TRUE vector is ignored.
#' @param collapse Logical (default FALSE). If TRUE the (m,p)-matrix will be
#' collapased into a (2*m+p)-vector containing the total number of matching
#' alleles.
#' @param wildcard Use the wildcard comparing.
#' @param wildcard.effect Compare result of wildcard and no wildcard.
#' @param wildcard.impose Force homozygouse profiles (aa) to have wildcard
#' (aF).
#' @param Rallele Implementation of 'Rare allele'designation matching.
#' @param threads The number of threads to use for performing comparisons in
#' parallel for increased computation time. Use 0 for using the same number as
#' the computer has CPU cores. NOTE: Only available on Linux and MacOS
#' operating systems.
#' @return Returns a matrix with the number of pairs
#' mathcing/partially-matching at (i,j)-loci.
#' @author James Curran and Torben Tvedebrink. The multicore/CPU implementation
#' was provided by Mikkel Meyer Andersen.
#' @examples
#' 
#'   \dontrun{
#'   data(dbExample)
#'   dbCompare(dbExample,hit=5,trace=TRUE)
#'   }
#' 
#' @importFrom utils head
#' 
#' @export dbCompare
dbCompare <- function(x, profiles = NULL, hit = 7, trace = TRUE, vector = FALSE, collapse = FALSE, 
  wildcard = FALSE, wildcard.effect = FALSE, wildcard.impose = FALSE, Rallele = FALSE, threads = 2) {
  
  if (all(wildcard, wildcard.effect)) 
    stop("Only one wildcard option is allowed to be 'TRUE' at any time")
  
  if (threads < 1) {
    stop("Must have at least one thread")
  }
  
  x.cp <- x  ## Keep copy of input
  
  isAmelogenin = function(x) {
    ## Looks after amelogenin loci - and drop these
    amel.candidate <- unlist(lapply(utils::head(x, n = 20), function(y) {
      any(grepl("^(x|y|X|Y|xx|xy|XX|XY)$", paste(y)))
    }))
    
    return(amel.candidate)
  }
  
  x <- x[, !isAmelogenin(x)]
  
  ## 
  stripID = function(x) {
    id = nid = NULL
    
    if (ncol(x)%%2 != 0) {
      ## Unequal number of columns - assumes first col is identifier
      id <- x[, 1, drop = TRUE]
      nid <- names(x)[1]
      x <- x[, -1]
    } else {
      id = 1:nrow(x)
      nid = NULL
    }
    
    return(list(x = x, id = id, nid = nid))
    
  }
  
  r = stripID(x)
  x = r$x
  id = r$id
  nid = r$nid
  
  ## number of loci
  numLoci <- ncol(x)/2
  
  toInteger = function(x) {
    ## Converts all alleles to a*10, e.g. 9 -> 90 and 9.3 -> 93 (to deal with .1, .2 and .3
    ## alleles)
    for (i in 1:ncol(x)) {
      if (inherits(x[[i]], "factor"))
        x[[i]] <- paste(x[[i]])
      x[[i]] <- as.numeric(x[[i]]) * 10
      class(x[[i]]) <- "integer"
    }
    
    return(x)
  }
  
  x = toInteger(x)
  
  ## Specific profile(s) that should be compared to profiles in x
  single <- 0
  
  if (!is.null(profiles)) {
    nx <- names(x)
    if (is.null(dim(profiles))) 
      dim(profiles) <- c(1, length(profiles))
    
    if (is.matrix(profiles)) 
      profiles <- as.data.frame(profiles)
    
    # amel.candidate <- unlist(lapply(head(profiles,n=min(nrow(profiles),20)),function(y)
    # any(grepl('^(x|y|X|Y|xx|xy|XX|XY)$',paste(y)))))
    profiles <- profiles[, !isAmelogenin(profiles)]
    
    if (ncol(profiles) == ncol(x)) {
      names(profiles) <- nx
    } else {
      id <- c(profiles[, 1, drop = TRUE], id)
      profiles <- structure(profiles[, -1], .Names = nx)
    }
    
    # for(i in 1:ncol(profiles)){ if(class(profiles[[i]])=='factor') profiles[[i]] <-
    # paste(profiles[[i]]) profiles[[i]] <- as.numeric(profiles[[i]])*10 class(profiles[[i]]) <-
    # 'integer' }
    profiles = toInteger(profiles)
    x <- rbind(profiles, x)
    single <- nrow(profiles)  ## Number of profiles to compare with x // input for C++ code
  }
  
  orderAlleles = function(x) {
    for (i in ((0:(numLoci - 1)) * 2 + 1)) {
      ## Looks for situations where A[1]>A[2]: C++ code assumes A[1]<=A[2]
      if (any(swap.index <- x[[i]] > x[[i + 1]])) {
        if (Rallele) 
          swap.index[x[[i]] == 990] <- FALSE  ## If R-allele (coded as 99*10, cf above) don't swap (order matters)
        x[swap.index, c(i, i + 1)] <- x[swap.index, c(i + 1, i)]
      }
      ## Force homs to have wildcards
      if (wildcard.impose) 
        x[[i]][x[[i]] == x[[i + 1]]] <- 0
    }
    return(x)
  }
  
  x = orderAlleles(x)
  
  if (!wildcard) {
    ## If no wildcard is allowed, but the database contains '0' which is equivalent of 'F' then
    ## these profiles are removed before comparison
    x0 <- x[(xnull <- apply(x[, -1], 1, function(y) any(y == 0))), ]
    x <- x[!xnull, ]
  }
  
  x.cp1 <- x  ## Keep copy of x
  res <- NULL
  if (threads > 1) {
    x <- do.call("paste", c(x = x, sep = "\t"))  ## Converts every line in the DB to a string separated by '\t'
    
    RcppParallel::setThreadOptions(numThreads = threads, stackSize = "auto")
    
    res <- compare_threaded(DB = x, numLoci = numLoci, bigHit = hit, trace = trace, single = single, 
                            useWildcard = wildcard, useWildcardEffect = wildcard.effect, 
                            useRallele = Rallele)
  } else {
    x <- do.call("paste", c(x = x, sep = "\t"))  ## Converts every line in the DB to a string separated by '\t'
    res <- compare(DB = x, numLoci = numLoci, bigHit = hit, trace = trace, single = single, 
      useWildcard = wildcard, useWildcardEffect = wildcard.effect, useRallele = Rallele)
  }
  
  stopifnot(!is.null(res))
  # print(str(res, 1))
  
  if (wildcard.effect) {
    result <- list(m = matrix(res$M, 2 * numLoci + 1, 2 * numLoci + 1, byrow = TRUE, dimnames = list(Genuine = 0:(2 * 
      numLoci), Wildcard = 0:(2 * numLoci))))
  } else {
    result <- list(m = matrix(res$M, numLoci + 1, numLoci + 1, byrow = TRUE, dimnames = list(match = 0:numLoci, 
      partial = 0:numLoci)))
  }
  
  call <- list(loci = numLoci, single = single, collapse = collapse, vector = vector, wildcard = c(wildcard, 
    wildcard.effect))
  if (length(res$row1) > 0) {
    result <- c(result, list(hits = data.frame(id1 = id[res$row1], id2 = id[res$row2], match = res$matches, 
      partial = res$partial, stringsAsFactors = FALSE)))
    if (wildcard) 
      result$hits <- cbind(result$hits, Fmatch = res$fmatches, Fpartial = res$fpartial)
    names(result$hits)[1:2] <- paste(nid, 1:2, sep = "")
    call <- c(call, list(hit = hit))
  } else result$hits <- NULL
  if (collapse) 
    result$m <- structure(dbCollapse(result$m), .Names = 0:(2 * numLoci)) else if (vector) 
    result$m <- structure(t(result$m)[up.tri(result$m)], .Names = t(outer(dimnames(result$m)[[1]], 
      dimnames(result$m)[[2]], paste, sep = "/"))[up.tri(result$m)])
  if (!wildcard) {
    if (nrow(x0) > 0) 
      result$excludedProfiles <- x0 else result$excludedProfiles <- "none"
  }
  attributes(result)$call <- call
  attributes(result)$class <- "dbcompare"
  result
}



#' Prints the summary matrix
#' 
#' Prints the summary matrix and possible 'big hits'.
#' 
#' Prints the summary matrix
#' 
#' @param x Summary matrix returned from dbcompare
#' @param ... ...
#' @return Prints the summary matrix and data frame with 'big hits'
#' @author James Curran and Torben Tvedebrink
#' @seealso dbCompare,plot.dbcompare
#' @examples
#' 
#'   \dontrun{
#'   data(dbExample)
#'   M = dbCompare(dbExample,hit=5)
#'   M
#'   }
#' 
#' @rawNamespace S3method(print, dbcompare)
#' @export print.dbcompare
print.dbcompare <- function(x, ...) {
  if (is.matrix(x$m)) {
    if (attributes(x)$call$wildcard[2]) 
      x$m[lower.tri(x$m)] <- NA else x$m[!up.tri(x$m)] <- NA
  }
  if (length(attributes(x)$names) == 1) {
    if (is.matrix(x$m)) 
      print.table(x$m, ...) else print.default(x$m, ...)
  } else {
    if (is.matrix(x$m)) {
      cat("Summary matrix\n")
      print.table(x$m, ...)
    } else {
      cat("Summary vector\n")
      print.default(x$m, ...)
    }
    if (!is.null(x$hits)) {
      cat(paste("\nProfiles with at least", attributes(x)$call$hit, "matching loci\n"))
      x$hits <- x$hits[with(x$hits, order(match, partial, decreasing = TRUE)), ]
      if ("Fmatch" %in% names(x$hits)) 
        x$hits <- x$hits[with(x$hits, order(match, partial, match, partial, decreasing = TRUE)), 
          ]
      rownames(x$hits) <- 1:nrow(x$hits)
      print.data.frame(x$hits, ...)
    }
  }
}

