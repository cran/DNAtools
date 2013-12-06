\name{pContrib}
\alias{pContrib}
\title{Compute the posterior probabilities for P(m|n0) for a given prior
  P(m) and observed vector n0 of locus counts} 
\description{
  Compute the posterior distribution \eqn{P(m|n_0)}{P(m|n0)}
  where \eqn{m}{m} ranges from 1 to \eqn{m_{max}}{\code{m.max}} and
  \eqn{n_0}{n0} is the observed locus counts.
}
\usage{
pContrib(n0, probs=NULL, m.prior=rep(1/m.max,m.max), m.max=8,
         pnoa=NULL, theta=0)
}
\arguments{
  \item{n0}{Vector of observed allele counts - same length as the number
    of loci.}
  \item{probs}{List of vectors with allele probabilities for each locus}
  \item{m.prior}{A vector with prior probabilities (summing to 1), where
    the length of \code{m.prior} determines the  plausible range of m}
  \item{m.max}{Derived from the length of \code{m.prior}, and if
    \code{m.prior=NULL} a uniform prior is speficied by \code{m.max}:
    \code{m.prior = rep(1/m.max,m.max)}.}
  \item{pnoa}{A named list of locus specific probabilities. Output from
    \code{pNoA} with \code{locuswise=TRUE}.}
  \item{theta}{The coancestery coefficient}
}
\details{
  Computes a vector \eqn{P(m|n0)}{P(m|n0)} evaluated over the plausible
  range \eqn{1,\dots,m_{max}}{1,...,m.max}.
}
\references{
  T Tvedebrink (2013). 'On the exact distribution of the number of
  alleles in DNA mixtures', International Journal of Legal Medicine:
  (under review). 
}
\value{
  Returns a vector P(m|n0) for m=1,...,m.max
}

\author{Torben Tvedebrink}

\examples{
  ## Simulate some allele frequencies:
  freqs <-  structure(replicate(10,
              { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE),.Names=paste("locus",1:10,sep="."))
  m <- 2
  n0 <- unlist(lapply(freqs,function(x)
         length(unique(sample(length(x),size=2*m,replace=TRUE,prob=x)))))
  ## Compute P(m|n0) for m=1,...,4 and the sampled n0
  pContrib(n0=n0,probs=freqs,m.max=4)
}

