\name{DNAtools-package}
\alias{DNAtools-package}
\alias{DNAtools}
\docType{package}
\title{
  Tools for analysing forensic genetic DNA databases
}
\description{
  Computational efficient tools for comparing all pairs of profiles in a
  DNA database. The expectation and covariance of the summary statistic
  is implemented for fast computing. Routines for estimating proportions
  of close related individuals are available. The use of wildcards (also
  called F-designation) is implemented. Dedicated functions ease
  plotting the results.
}
\details{
\tabular{ll}{
Package: \tab DNAtools\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-08-25\cr
License: \tab GPL (>= 2)\cr
}
dbCompare: Compares make all n(n-1)/2 pairwise comparisons between
profiles of a database with n DNA profiles.
dbExpect: Computes the expected number of matching and partial matching
loci for a given number of profiles in a database.
dbVariance: Calculates the associated covariance matrix.
}
\author{
Torben Tvedebrink <tvede@math.aau.dk> and James Curran
<j.curran@auckland.ac.nz> wrote the package. Mikkel Meyer Andersen
<mikl@math.aau.dk> assisted on the multicore/CPU implementation.

Maintainer: Torben Tvedebrink <tvede@math.aau.dk>
}
\references{
  Tvedebrink T, JM Curran, PS Eriksen, HS Mogensen and N Morling (2012).
  Analysis of matches and partial-matches in a Danish STR data set.
  Forensic Science International: Genetics, 6(3): 387-392.

  Read the vignette: \code{vigette("DNAtools")}
}
\keyword{ Forensic genetics }
\examples{
  \dontrun{
  data(dbExample)
  dbCompare(dbExample,hit=5,trace=TRUE)
  }
}
