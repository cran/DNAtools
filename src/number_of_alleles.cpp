#include <Rcpp.h>

#include "class_probsObj.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @export
// [[Rcpp::export]]
NumericVector convolve(NumericMatrix x) {
  int L = x.nrow();
  int m = x.ncol() / 2;
  NumericMatrix P(L, 2 * m * L);
  int i, j, l; 
  double q;
  
  for(j = 0; j < 2 * m; j++) 
    P(0, j) = x(0, j);
  
  for(l = 1; l < L; l++) { // Loop over loci
    /*
    end bound is for l = L-1, 
    then 2 * m * (l + 1) = 2 * m * (L-1 + 1) = 2*m*L
    */
    //for(j = 0 ; j < 2 * m * (l + 1); j++) { // Loop over min and max value
    for(j = 0 ; j < 2 * m * (l + 1) - 1; j++) { // Loop over min and max value
      q = 0.0;
      
      for(i = 0; i < 2 * m; i++){
        if(j - i >= 0) 
          q += P(l - 1, j - i) * x(l , i);
      }
      P(l, j + 1) = q; // FIXME/ERROR: mikl: j + 1; upper bound changed from 2 * m * (l + 1) to 2 * m * (l + 1) - 1
    }
  }
  
  NumericVector result = P(L - 1, _);
  
  Rcpp::CharacterVector nms(result.size());
  
  for(int i = 0; i < result.size(); i++){
    nms(i) = std::to_string(i + 1);
  }
	result.attr("names") = nms;
  
  return result;
}

//' @export
// [[Rcpp::export]]
NumericVector Pnm_locus(int m, double theta, NumericVector alleleProbs){
  probsObj P(alleleProbs, theta);
  NumericVector result = P.calcProbs(m);

  //Rprintf("%s\n", P.printLookup().c_str());
  return result;
}

// [[Rcpp::export]]
NumericMatrix Pnm_all_cpp(int numContrib, double theta, List loci){
  probsObj P(theta);
  NumericMatrix result = P.calcProbs(numContrib, loci);  
  return result;
}

// [[Rcpp::export]]
List generateCompositions(int numContributors) {
  probsObj P;
  return P.getCompositions(numContributors);
}



