\name{plot.dbOptim}
\alias{plot.dbOptim}
\alias{points.dbOptim}
\alias{lines.dbOptim}
\title{Plots the fitted object function for estimated familial
  relationships in the database and theta.} 
\description{
  Plots the minimised object function for included values of theta
}
\usage{
\method{plot}{dbOptim}(x,type="l",...)
\method{points}{dbOptim}(x,type="p",...)
\method{lines}{dbOptim}(x,type="l",...)
}
\arguments{
  \item{x}{Object returned by optim.relatedness}
  \item{type}{The type of plot character ("l"=line, "p"=points, ...),
    see 'par' for more details}
  \item{...}{Other plot options}
}
\details{
  Plots the object function
}
\value{
  A plot of the object function
}

\seealso{optim.relatedness}

\author{James Curran and Torben Tvedebrink}

\examples{
  \dontrun{
  ## Simulate some allele frequencies:
  freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
              simplify=FALSE)
  ## Load the sample database:
  data(dbExample)
  obs <- dbCompare(dbExample,trace=FALSE)$m
  C3 <- optim.relatedness(obs,theta0=0.0,theta1=0.03,probs=freqs,
          objFunction="C3",max.bisect=30,trace=TRUE)
  plot(C3)
  }
}

