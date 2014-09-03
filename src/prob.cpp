#include <Rcpp.h>
#include "DNTRare.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector Prob(CharacterVector vstrCombs, NumericVector q, NumericVector R, 
                  double r, double t){
 
  DNTRare dnt(q, R, r, t);
  NumericVector vResult = dnt.prob(as<vector<string> >(vstrCombs));
  
  return vResult;
}
