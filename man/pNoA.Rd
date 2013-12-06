\name{pNoA}
\alias{pNoA}
\alias{p.numberofalleles}
\title{The exact distribution of the number of alleles in a m-person DNA mixture}
\description{
  Computes the exact distribution of the number of alleles in a m-person
  DNA mixture typed with STR loci. For a m-person DNA mixture it is
  possible to observe 1,...,2*m*L alleles, where L is the total number
  of typed STR loci. The method allows incorporation of the
  subpopulation correction, the so-called
  \eqn{\theta}{theta}-correction, to adjust for shared ancestry. If
  needed, the locus-specific probabilities can be obtained using 
  the \code{locuswise} argument.
}
\usage{
pNoA(probs, m=2, theta = 0, noa.tabs=NULL, locuswise=FALSE)
}
\arguments{
  \item{probs}{List of vectors with allele probabilities for each locus}
  \item{m}{The number of contributors}
  \item{theta}{The coancestery coefficient}
  \item{noa.tabs}{If \code{noa.tabs} has been computed by the
    \code{noaTabs}-function for \code{m}, then this will speed up the
    computations.} 
  \item{locuswise}{Logical. If TRUE the locuswise probabilities will be
    returned. Otherwise, the probability over all loci is returned.}
}
\details{
  Computes the exact distribution of the number of alleles for a
  m-person DNA mixture.
}
\references{
  T Tvedebrink (2013). 'On the exact distribution of the number of
  alleles in DNA mixtures', International Journal of Legal Medicine:
  (under review). 
}
\value{
  Returns a vector of probabilities, or a list of locuswise probability vectors.
}

\author{Torben Tvedebrink}

\examples{
  ## Simulate some allele frequencies:
  freqs <-  structure(replicate(10,
              { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE),.Names=paste("locus",1:10,sep="."))
  ## Compute alpha_3
  noa.tab3 <- noaTabs(m=3)

  ## Compute P(N(m=3)=n), n=1,...,2*L*m, where L=10 here
  pNoA(freqs,m=3,theta=0,noa.tabs=noa.tab3)
  ## Same, but locuswise results
  pNoA(freqs,m=3,theta=0,noa.tabs=noa.tab3,locuswise=TRUE)
}

