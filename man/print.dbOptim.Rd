\name{print.dbOptim}
\alias{print.dbOptim}
\title{Prints the results from optim.relatedness()}
\description{
  Prints the evaluated functions for the object function, best estimate
  of alpha and possibly list of variances. 
}
\usage{
print.dbOptim(x,var.list=FALSE,...)
}
\arguments{
  \item{x}{Object returned by optim.relatedness()}
  \item{var.list}{Logical. Whether the (long) list of variance
    components should be printed to the screen.}
  \item{...}{...}
}
\details{
  Prints the summary details of the fit
}
\value{
  A dataframe with [theta,value] and a vector of fitted alpha parameters 
}

\seealso{optim.relatedness}

\author{James Curran and Torben Tvedebrink}

\examples{
  \dontrun{
  ## Simulate some allele frequencies:
  freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)}, simplify=FALSE)
  ## Load the sample database:
  data(dbExample)
  obs <- dbCompare(dbExample,trace=FALSE)$m
  C3 <- optim.relatedness(obs,theta0=0.0,theta1=0.03,probs=freqs,objFunction="C3",max.bisect=30,trace=TRUE)
  print(C3)
  }
}
