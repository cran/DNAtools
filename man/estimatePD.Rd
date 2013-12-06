\name{estimatePD}
\alias{estimatePD}
\title{Estimate the drop-out probability based on number of alleles}
\description{
  An inferior may to estimate the drop-out probability compared to using
  the peak heights from the electropherogram. However, to compare the
  performance with Gill et al. (2007) this implements a theoretical
  approach based on their line of arguments.
}
\usage{
estimatePD(n0,m,pnoa=NULL,probs=NULL,theta=0,noa.tabs=NULL,
           locuswise=FALSE)
}
\arguments{
  \item{n0}{Vector of observed allele counts - same length as the number
    of loci.}
  \item{m}{The number of contributors}  
  \item{pnoa}{The vector of P(N(m)=n) for n=1,...,2Lm, where L is the
    number of loci and m is the number of contributors.}
  \item{probs}{List of vectors with allele probabilities for each locus}
  \item{theta}{The coancestery coefficient}
  \item{noa.tabs}{If \code{noa.tabs} has been computed by the
    \code{noaTabs}-function for \code{m}, then this will speed up the
    computations.}
  \item{locuswise}{Logical. Indicating whether computations should be
    done locuswise.}
}
\details{
  Computes the P(D) that maximises equation (10) in Tvedebrink (2013).
}
\references{
  Gill, P., A. Kirkham, and J. Curran (2007).
  LoComatioN: A software tool for the analysis of low copy number DNA profiles.
  Forensic Science International 166(2-3): 128 - 138.

  T Tvedebrink (2013). 'On the exact distribution of the number of
  alleles in DNA mixtures', International Journal of Legal Medicine:
  (under review). 
}
\value{
  Returns the MLE of P(D) based on equation (10) in Tvedebrink (2013)
}

\author{Torben Tvedebrink}

\examples{
  ## Simulate some allele frequencies:
  freqs <-  structure(replicate(10,
              { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE),.Names=paste("locus",1:10,sep="."))
  ## Assume 15 alleles are observed in a 2-person DNA mixture with 10 loci:
  estimatePD(n0=15,m=2,probs=freqs)
}

