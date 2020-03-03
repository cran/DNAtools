#ifndef _COMPARE_UTILS_H
#define _COMPARE_UTILS_H 1

using namespace std;

#include <Rcpp.h>

#include "profile.cpp"

vector<Profile *> readProfiles(const Rcpp::StringVector &DB, int nProfiles, int nLoci);

Rcpp::List prepReturnList(Rcpp::IntegerVector &m, vector<int>& vnRow1, vector<int>& vnRow2, vector<int>& vnMatch, 
                          vector<int>& vnPartial, vector<int>& vnFmatch, vector<int>& vnFpartial);

#endif
