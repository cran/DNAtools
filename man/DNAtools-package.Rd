\name{DNAtools-package}
\alias{DNAtools-package}
\alias{DNAtools}
\docType{package}
\title{
Tools for analysing forensic genetic DNA data
}
\description{
Functions for making a database comparison exercise where each DNA
profile is compared to any other in the database. Also has functions for
optimising parameters of the population genetic model and for simulating
DNA databases. 
}
\details{

Index:
\preformatted{
dbCompare                                                   Compare DNA profiles
dbExample                                                   Simulated database with 1,000 individuals
dbExpect                                                    Expected value of cell counts in DNA database comparison
dbSimulate                                                  Simulate a DNA database
dbVariance                                                  Covariance matrix of cell counts in DNA database comparison
optim.relatedness                                           Estimate theta and the fraction of comparisons between close relatives
freqEst                                                     Simple allele frequency estimation
plot.dbOptim                                                Plots the fitted object function for estimated familial relationships
plot.dbcompare                                              Plots the summary matrix
print.dbOptim                                               Prints the results from optim.relatedness()
print.dbcompare                                             Prints the summary matrix
}
}
\author{
Torben Tvedebrink and James Curran

Maintainer: Torben Tvedebrink <tvede@math.aau.dk>
}
\references{
  T Tvedebrink, PS Eriksen, J Curran, HS Mogensen, N Morling. 'Analysis
  of matches and partial-matches in Danish DNA reference profile
  database'. Forensic Science International: Genetics, 2011.
}
\keyword{ package }
\examples{
data(dbExample)
freq <- freqEst(dbExample)
M <- dbCompare(dbExample,trace=FALSE)
plot(M)
M.expect <- dbExpect(probs=freq,theta=0.0,n=nrow(dbExample),vector=TRUE)
points(M.expect,col=2,pch=4)
legend("topright",bty="n",pch=c(1,4),col=c(1,2),c("Observed", "Expected"))
}
