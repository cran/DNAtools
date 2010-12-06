\name{freqEst}
\alias{freqEst}
\title{Simple allele frequency estimation}
\description{
  Estimates allele frequencies from a database with DNA profiles
}
\usage{
freqEst(x)
}
\arguments{
  \item{x}{A database of the form ["id","locus1 allele1","locus1
    allele2",...,"locusN allele 1","locusN allele2"].}
}
\details{
  Computes the allele frequencies for a given database.
}
\value{
  Returns a list of probability vectors - one vector for each locus.
}
  
\author{James Curran and Torben Tvedebrink}

\examples{
  data(dbExample)
  freqEst(dbExample)
}


