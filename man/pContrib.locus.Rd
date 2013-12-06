\name{pContrib.locus}
\alias{pContrib.locus}
\title{Compute the posterior probabilities for P(m|n0) for a given prior P(m)}
\description{
  Compute a matrix of posterior probabilties \eqn{P(m|n_0)}{P(m|n0)}
  where m ranges from 1 to \eqn{m_{max}}{\code{m.max}}, and \eqn{n_0}{n0}
  is \eqn{0,\dots,2m_{max}}{0,..,\code{m.max}}. This is done by
  evaluating \eqn{P(m|n_0)=P(n_0|m)P(m)/P(n)}{P(m|n0)=P(n0|m)P(m)/P(n)},
  where \eqn{P(n_0|m)}{P(n0|m)} is evaluated by \code{pNoA}. 
}
\usage{
pContrib.locus(prob=NULL, m.prior=NULL, m.max=8, pnoa.locus=NULL,
               theta=0)
}
\arguments{
  \item{prob}{Vectors with allele probabilities for the specific locus}
  \item{m.prior}{A vector with prior probabilities (summing to 1), where
    the length of \code{m.prior} determines the  plausible range of m}
  \item{m.max}{Derived from the length of \code{m.prior}, and if
    \code{m.prior=NULL} a uniform prior is speficied by \code{m.max}:
    \code{m.prior = rep(1/m.max,m.max)}.}
  \item{pnoa.locus}{A named vector of locus specific probabilities
    \eqn{P(N(m)=n), n=1,\dots,2m}{P(N(m)=n), n=1,..,2m}.}
  \item{theta}{The coancestery coefficient}
}
\details{
  Computes a matric of \eqn{P(m|n0)}{P(m|n0)} values for a specific locus.
}
\references{
  T Tvedebrink (2013). 'On the exact distribution of the number of
  alleles in DNA mixtures', International Journal of Legal Medicine. 
}
\value{
  Returns a matrix [P(m|n0)] for m=1,...,m.max and n0=1,...,2m.max.
}

\author{Torben Tvedebrink}

\examples{
  ## Simulate some allele frequencies:
  freqs <-  structure(replicate(10,
              { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE),.Names=paste("locus",1:10,sep="."))
  ## Compute P(m|n0) for m=1,...,5 and n0=1,...10 for the first locus:
  pContrib.locus(prob=freqs[[1]],m.max=5)
}

