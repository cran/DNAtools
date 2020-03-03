#include <vector>
#include <string>
#include <Rcpp.h>
#include "class_DNTRare.cpp"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Prob(CharacterVector vstrCombs, NumericVector q, NumericVector R, 
                  double r, double t){
 
  DNTRare dnt(q, R, r, t);
  NumericVector vResult = dnt.prob(as<vector<string> >(vstrCombs));
  
  return vResult;
}
