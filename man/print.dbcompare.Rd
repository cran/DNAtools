\name{print.dbcompare}
\alias{print.dbcompare}
\title{Prints the summary matrix}
\description{
  Prints the summary matrix and possible "big hits". 
}
\usage{
\method{print}{dbcompare}(x,...)
}
\arguments{
  \item{x}{Summary matrix returned from dbcompare}
  \item{...}{...}
}
\details{
  Prints the summary matrix
}
\value{
  Prints the summary matrix and data frame with "big hits"
}

\seealso{\code{\link{dbCompare}}, \code{\link{plot.dbcompare}}}

\author{James Curran and Torben Tvedebrink}

\examples{
  \dontrun{
  data(dbExample)
  M = dbCompare(dbExample,hit=5)
  M
  }
}


