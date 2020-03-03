#' Expected value of cell counts in DNA database comparison
#' 
#' Computes the expected number of cell counts when comparing DNA profiles in a
#' DNA database. For every pair of DNA profiles in a database the number of
#' matching and partial matching loci is recorded. A match is declared if the
#' two DNA profiles coincide for both alleles in a locus and a partial-match is
#' recorded if only one allele is shared between the profiles. With a total of
#' L loci the number of matching loci is 0,...,L and partial number of matches
#' is 0,...,L-m, where m is the number of matching loci.
#' 
#' Computes the expected cell counts using a recursion formula. See Tvedebrink
#' et al (2011) for details.
#' 
#' @param probs List of vectors with allele probabilities for each locus
#' @param theta The coancestery coefficient
#' @param k The vector of identical-by-descent probabilities, k=(k2,k1,k0),
#' where for full-siblings k=c(1,2,1)/4. The default is k=c(0,0,1) refering to
#' unrelated individuals.
#' @param n Number of DNA profiles in the database
#' @param r The probability assigned to the rare alleles (see rare allele
#' matching). If a vector must be of same length as \code{probs}.
#' @param R The probability assigned to alleles shorter or longer than allelic
#' ladder (see rare allele matching). If a vector must be of length 1 or 2, and
#' if a list it must be same length as \code{probs}.
#' @param round Whether or not the results should be rounded or not
#' @param na Whether or not the off-elements should be returned as 0 or NA
#' @param vector Whether or not the result should be returned as a matrix or
#' vector. Note if 'collapse' is TRUE vector is ignored.
#' @param collapse Logical (default FALSE). If TRUE the (m,p)-matrix will be
#' collapased into a (2*m+p)-vector containing the total number of matching
#' alleles.
#' @param wildcard Should wildcards be used?
#' @param no.wildcard Should 'w' wildcards be used?
#' @param rare.allele Should rare allele matching be used?
#' @param no.rare.allele Should 'r' rare allele loci be used?
#' @return Returns a matrix (or vector, see above) of expected cell counts.
#' @author James Curran and Torben Tvedebrink
#' @references T Tvedebrink, PS Eriksen, J Curran, HS Mogensen, N Morling.
#' 'Analysis of matches and partial-matches in Danish DNA reference profile
#' database'. Forensic Science International: Genetics, 2011.
#' @examples
#' 
#'   \dontrun{
#'   ## Simulate some allele frequencies:
#'   freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
#'               simplify=FALSE)
#'   ## Compute the expected number for a DB with 10000 profiles:
#'   dbExpect(freqs,theta=0,n=10000)
#'   } 
#' 
#' @export dbExpect
dbExpect <- function(probs, theta = 0, k = c(0, 0, 1), n = 1, r = 0, R = 0, round = FALSE, na = TRUE, 
  vector = FALSE, collapse = FALSE, wildcard = FALSE, no.wildcard = NULL, rare.allele = FALSE, 
  no.rare.allele = NULL) {
  if (length(theta) > 1) {
    return(lapply(theta, function(t) dbExpect(probs = probs, theta = t, k = k, n = n, r = r, 
      R = R, round = round, na = na, vector = vector, collapse = collapse, wildcard = wildcard, 
      no.wildcard = no.wildcard, rare.allele = rare.allele, no.rare.allele = no.rare.allele)))
  }
  
  if (rare.allele & wildcard) {
    stop("Only one of 'wildcard' and 'rare.allele' can be TRUE.\nOtherwise set 'no.wildcard' and 'no.rare.allele'.")
  }
  
  if (wildcard) {
    f <- "F"
  } else if (rare.allele) {
    if (any(c(unlist(r), unlist(R)) > 0)) {
      f <- "R"  ## If one or more thresholds or tail probabilities are set
    } else {
      f <- ""  ## If no probability is assigned the rare alleles
    }
  } else {
    f <- ""
  }
  
  if (!is.list(probs) && is.vector(probs)) {
    return(get(paste(f, "Ps", sep = ""))(probs, t = theta, k = k, r = r, R = R))  ## Previous line handles if more thetas are provided 
  }
  
  ## probs is a list of vectors with each vector being the allele probabilities for a given
  ## locus
  S <- length(probs)
  if (any(unlist(r) > 0)) {
    if (!(length(r) == 1 | length(r) == S)) {
      stop("Length of 'r' must be equal to the number of loci.")
    }
    
    if (length(r) == 1) {
      r <- rep(r, S)
    }
  }
  
  if (any(unlist(R) > 0)) {
    if (!(length(unlist(R)) %in% c(1, 2, S, 2 * S))) {
      stop("Length of 'R' not correct.")
    }
    
    if (length(R) == 1 && is.list(R)) {
      R <- replicate(S, R, simplify = FALSE)
    } else if (!is.list(R)) {
      R <- replicate(S, R, simplify = FALSE)
    }
  }
  
  if (is.character(k)) {
    if (!is.na(match(toupper(k), c("UN", "AV", "FS", "FC", "PC")))) {
      k <- list(UN = c(0, 0, 1), FC = c(0, 1, 3)/4, AV = c(0, 1, 1)/2, PC = c(0, 1, 0), 
        FS = c(1, 2, 1)/4)[[toupper(k)]]
    } else {
      stop(paste(paste("The value of 'k' (", k, ") is not defined - has to be vector c(k2,k1,k0) or string.", 
        sep = ""), "Options are: 'UN' (unrelated), 'FS' (full-siblings), 'AV' (avuncular), 'FC' (first cousins) or 'PC' (parent-child)", 
        sep = "\n"))
    }
  }
  
  if (length(r) > 1) {
    if (length(R) != length(r)) {
      R <- replicate(length(r), R, simplify = FALSE)
    }
    p <- lapply(as.list(1:S), function(z) {
      RPs(probs[[z]], r = r[z], t = theta, k = k, R = R[[z]])
    })
  } else {
    p <- lapply(probs, get(paste(f, "Ps", sep = "")), t = theta, k = k, r = r, R = R)
  }
  
  q <- do.call("cbind", p)
  res <- rare(q)
  ## In case of mixed matching
  if (!is.null(no.wildcard) | !is.null(no.rare.allele)) {
    if (is.null(no.wildcard)) {
      no.wildcard <- 0
    }
    
    if (is.null(no.rare.allele)) {
      no.rare.allele <- 0
    }
    
    if ((no.wildcard + no.rare.allele) > S) {
      stop("Too many wildcard and rare allele loci compared to total number of loci")
    }
    
    if (no.wildcard == 0 & no.rare.allele == 0) {
      return(dbExpect(probs = probs, theta = theta, k = k, n = n, r = r, R = R, round = round, 
        na = na, vector = vector, collapse = collapse, wildcard = FALSE, rare.allele = FALSE))
    } else if (no.wildcard == S || no.rare.allele == S) {
      ## If w+r equals L // No exact matching If all wildcard -> no rare
      if (no.wildcard == S) {
        return(dbExpect(probs = probs, theta = theta, k = k, n = n, r = r, R = R, round = round, 
          na = na, vector = vector, collapse = collapse, wildcard = TRUE, rare.allele = FALSE))
      } else {
        ## 
        return(dbExpect(probs = probs, theta = theta, k = k, n = n, r = r, R = R, round = round, 
          na = na, vector = vector, collapse = collapse, wildcard = FALSE, rare.allele = TRUE))
      }
    } else {
      ## If mixture
      M <- permAll(rep(0:2, c(no.wildcard, no.rare.allele, S - (no.wildcard + no.rare.allele))))
      nM <- nrow(M)
      PPs <- do.call("cbind", lapply(1:length(probs), function(x) Ps(p = probs[[x]], t = theta, 
        k = k, r = r[x], R = R[[x]])))
      FFs <- do.call("cbind", lapply(1:length(probs), function(x) FPs(p = probs[[x]], 
        t = theta, k = k, r = r[x], R = R[[x]])))
      RRs <- do.call("cbind", lapply(1:length(probs), function(x) RPs(p = probs[[x]], 
        t = theta, k = k, r = r[x], R = R[[x]])))
      res <- rare(cbind(FFs[, M[1, ] == 0], RRs[, M[1, ] == 1], PPs[, M[1, ] == 2]))/nM
      for (i in 2:nM) {
        res <- res + rare(cbind(FFs[, M[i, ] == 0], RRs[, M[i, ] == 1], PPs[, M[i, ] == 
          2]))/nM
      }
    }
  }
  
  if (n > 1) {
    N <- choose(n, 2)
  } else {
    N <- 1
  }
  
  if (round) {
    res <- round(res * N)
  } else {
    res <- res * N
  }
  
  if (na) {
    res[!up.tri(res)] <- NA
  }
  
  if (collapse) {
    return(dbCollapse(res))
  }
  
  if (vector) {
    res <- t(res)[up.tri(res)]
    names(res) <- dbCats(S, vector = TRUE)
  } else {
    dimnames(res) <- list(match = 0:S, partial = 0:S)
  }
  res
}

