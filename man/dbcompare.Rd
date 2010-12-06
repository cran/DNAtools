\name{dbCompare}
\alias{dbCompare}
\title{Compare DNA profiles}
\description{
  Compare DNA profiles
}
\usage{
dbCompare(x, profile=NULL, hit=7, trace=TRUE, vector=FALSE)
}
\arguments{
  \item{x}{Database with DNA profiles. The database format is expected
    to be a data frame with each column containing an allelic number
    such that for each DNA marker there are two columns in the data
    frame. See \code{data(dbExample)} for an example of the format.}
  \item{profile}{A specific profile to be compared with all profiles in
    the database. Input is a vector of same length as a row in the
    database \code{x}.}
  \item{hit}{The number of matching loci for further investigation}
  \item{trace}{Shows a progress bar}
  \item{vector}{Logical. Whether the result should be returned as vector
    or a matrix}
}
\details{
  Computes the distance between DNA profiles in terms of matching and
  partially-matching STR loci.
}
\value{
  Returns a matrix with the number of pairs mathcing/partially-matching at (i,j)-loci.
}
  
\author{James Curran and Torben Tvedebrink}

\examples{
  \dontrun{
  data(dbExample)
  dbCompare(dbExample,hit=5,trace=TRUE)
  }
}


