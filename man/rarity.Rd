\name{dbExpect}
\alias{dbExpect}
\title{Expected value of cell counts in DNA database comparison}
\description{
  Computes the expected number of cell counts when comparing DNA
  profiles in a DNA database. For every pair of DNA profiles in a
  database the number of matching and partial matching loci is
  recorded. A match is declared if the two DNA profiles coincide for
  both alleles in a locus and a partial-match is recorded if only one
  allele is shared between the profiles. With a total of L loci the
  number of matching loci is 0,...,L and partial number of matches is
  0,...,L-m, where m is the number of matching loci.
}
\usage{
dbExpect(probs,theta=0.03,k=c(0,0,1),n=1,round=FALSE,na=TRUE,vector=FALSE,collapse=FALSE)
}
\arguments{
  \item{probs}{List of vectors with allele probabilities for each locus}
  \item{theta}{The coancestery coefficient}
  \item{k}{The vector of identical-by-descent probabilities,
    k=(k2,k1,k0), where for full-siblings k=c(1,2,1)/4. The default is
    k=c(0,0,1) refering to unrelated individuals.}
  \item{n}{Number of DNA profiles in the database}
  \item{round}{Whether or not the results should be rounded or not}
  \item{na}{Whether or not the off-elements should be returned as 0 or NA}
  \item{vector}{Whether or not the result should be returned as a matrix
    or vector. Note if 'collapse' is TRUE vector is ignored.}
  \item{collapse}{Logical (default FALSE). If TRUE the
    (m,p)-matrix will be collapased into a (2*m+p)-vector containing the
    total number of matching alleles.}
}
\details{
  Computes the expected cell counts using a recursion formula. See
  Tvedebrink et al (2011) for details.}
\references{
  T Tvedebrink, PS Eriksen, J Curran, HS Mogensen, N Morling. 'Analysis
  of matches and partial-matches in Danish DNA reference profile
  database'. Forensic Science International: Genetics, 2011.
}
\value{
  Returns a matrix (or vector, see above) of expected cell counts.
}

\author{James Curran and Torben Tvedebrink}

\examples{
  \dontrun{
  ## Simulate some allele frequencies:
  freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)}, simplify=FALSE)
  ## Compute the expected number for a DB with 10000 profiles:
  dbExpect(freqs,theta=0.03,n=10000)
  } 
}

