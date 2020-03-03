#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>

using namespace std;

#include <Rcpp.h>

#include "compare-utils.h"

vector<Profile *> readProfiles(const Rcpp::StringVector &DB, int nProfiles, int nLoci){
  int t = 0;
  string strLine;
  vector<Profile *> vpProfiles;
  Profile *pProfile;
  
  while(t < nProfiles){
    strLine = DB(t);
    pProfile = new Profile(strLine, nLoci);
    vpProfiles.push_back(pProfile);
    t++;
  }
  
  return vpProfiles;
}

Rcpp::List prepReturnList(Rcpp::IntegerVector &m, vector<int>& vnRow1, vector<int>& vnRow2, vector<int>& vnMatch, 
                          vector<int>& vnPartial, vector<int>& vnFmatch, vector<int>& vnFpartial){
  
 Rcpp::List rl;
  
  //    Rprintf("vnRow1 %d\n",vnRow1.size());
  //    Rprintf("vnRow2 %d\n",vnRow2.size());
  //    Rprintf("vnMatch %d\n",vnMatch.size());
  //    Rprintf("vnPartial %d\n",vnPartial.size()); 
  //    Rprintf("vnFMatch %d\n",vnFmatch.size());
  //    Rprintf("vnFpartial %d\n",vnFpartial.size());
  
  int pnSizes[6];
  pnSizes[0] = (int)vnRow1.size();
  pnSizes[1] = (int)vnRow2.size();
  pnSizes[2] = (int)vnMatch.size();
  pnSizes[3] = (int)vnPartial.size();
  pnSizes[4] = (int)vnFmatch.size();
  pnSizes[5] = (int)vnFpartial.size();
  
  sort(pnSizes, pnSizes + 6);
  
  if(pnSizes[0] != pnSizes[5]){
    Rprintf("Warning: different result vector sizes in prepReturnList. This will cause problems\n");
    
    pnSizes[0] = (int)vnRow1.size();
    pnSizes[1] = (int)vnRow2.size();
    pnSizes[2] = (int)vnMatch.size();
    pnSizes[3] = (int)vnPartial.size();
    pnSizes[4] = (int)vnFmatch.size();
    pnSizes[5] = (int)vnFpartial.size();
    
    const char *szNames[] = {"vnRow1", "vnRow2", "vnMatch", "vnPartial", "vnFmatch", "vnFpartial"};
    for(int i = 0; i < 6; i++){
      Rprintf("%s %d\n", (char *)szNames[i], pnSizes[i]);
    }
  }
  
  
  int nMatchlength = (int)vnRow1.size();
  
  Rcpp::IntegerVector row1(vnRow1.size());
  Rcpp::IntegerVector row2(vnRow2.size());
  Rcpp::IntegerVector matches(vnMatch.size());
  Rcpp::IntegerVector partial(vnPartial.size());
  Rcpp::IntegerVector fmatches(vnFmatch.size());
  Rcpp::IntegerVector fpartial(vnFpartial.size());
  
  for(int i = 0; i < nMatchlength; i++){
    row1[i] = vnRow1[i]; 
    row2[i] = vnRow2[i]; 
    matches[i] = vnMatch[i]; 
    partial[i] = vnPartial[i]; 
    fmatches[i] = vnFmatch[i]; 
    fpartial[i] = vnFpartial[i]; 
  }
  
  rl["M"] = m;
  rl["row1"] = row1;
  rl["row2"] = row2;
  rl["matches"] = matches;
  rl["partial"] = partial;
  rl["fmatches"] = fmatches;
  rl["fpartial"] = fpartial;
  
  return rl;
}


// [[Rcpp::export]]
Rcpp::IntegerVector score_rcpp(const Rcpp::IntegerVector& prof1, 
                               const Rcpp::IntegerVector& prof2, 
                               int numLoci, 
                               bool useWildCard, 
                               bool useRareAllele){
  
  int *pnProf1 = &(Rcpp::as<std::vector<int> >(prof1)[0]); // this should replace INTEGER(prof1)
  int *pnProf2 = &(Rcpp::as<std::vector<int> >(prof2)[0]);
  
  int nLoci = numLoci;
  bool bWildCard = useWildCard;
  bool bRareAllele = useRareAllele;
  
  Profile *pProf1 = new Profile(pnProf1, nLoci);
  Profile *pProf2 = new Profile(pnProf2, nLoci);
  
  //Rprintf("Profile 1:\n%s\nProfile 2:\n%s\n", pProf1->toString().c_str(), pProf2->toString(true).c_str());
  
  vector<int> vnScore(nLoci);
  unsigned long m2, m1, m0, fm2, fm1;
  m2 = m1 = m0 = fm2  = fm1 = 0;
  
  pProf1->compare(pProf2, m2, m1, m0, fm2, fm1, bWildCard, bRareAllele, &vnScore);
  delete pProf1;
  delete pProf2;
  
  Rcpp::IntegerVector result;
  vector<int>::iterator i = vnScore.begin();
  while(i!=vnScore.end()){
    result.push_back(*i);
    i++;
  }
  return result;
}


