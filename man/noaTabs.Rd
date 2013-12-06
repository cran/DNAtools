\name{noaTabs}
\alias{noaTabs}
\title{Computes the alpha_m variants and their weights.}
\description{
  Computes the \eqn{\alpha_m}{alpha_m}-vectors and their weights,
  \eqn{c(\alpha_m)}{c(alpha_m)}, in order to compute the exact
  distribution of the number of alleles in a m-person DNA mixture.
}
\usage{
noaTabs(alpha=NULL,m=2,weight=1)
}
\arguments{
  \item{alpha}{If \code{noaTabs} has been evaluated for \code{n<m}, then
    providing this object to \code{noaTabs} reduces the number of
    recursions needed to evaluate for \code{m}.} 
  \item{m}{The number of contributors}
  \item{weight}{See the paper for details. Usually, this should be set
    to 1}
}
\details{
  Computes the \eqn{\alpha_m}{alpha_m}-vectors and their weights,
  \eqn{c(\alpha_m)}{c(alpha_m)}, by recursion over the number of
  contributors \eqn{m}{m}. That is, \eqn{\alpha_m}{alpha_m} is obtained
  based on \eqn{\alpha_{m-1}}{alpha_(m-1)} -- see paper for details.
}
\references{
  T Tvedebrink (2013). 'On the exact distribution of the number of
  alleles in DNA mixtures', International Journal of Legal Medicine:
  (under review). 
}
\value{
  Returns a named vector of counts, where the names denotes the
  \eqn{\alpha_m}{alpha_m}-vectors and the countes is the associated
  weights, \eqn{c(\alpha_m)}{c(alpha_m)}.
}

\author{Torben Tvedebrink}

\examples{
  ## Simulate some allele frequencies:
  freqs <-  structure(replicate(10,
              { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE),.Names=paste("locus",1:10,sep="."))
  ## Compute alpha_2
  noa.tab2 <- noaTabs(m=2)
  ## Use alpha_2 to compute alpha_3
  noa.tab3 <- noaTabs(alpha=noa.tab2, m=3)

  ## Compute P(N(m=3)=n), n=1,...,2*L*m, where L=10 here
  pNoA(freqs,m=3,theta=0,noa.tabs=noa.tab3)
  ## Same, but locuswise results
  pNoA(freqs,m=3,theta=0,noa.tabs=noa.tab3,locuswise=TRUE)
}

