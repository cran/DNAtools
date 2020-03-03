/*########################################################################*/
/* File: compare.cpp                                                     */
/* Authors: James M. Curran, Torben Tvedebrink, and                       */
/*          Mikkel M. Andersen                                            */
/*                                                                        */
/* Version history                                                        */
/* ---------------------------------------------------------------------  */
/* Version  Date        Changes                          Author           */
/* -------  ----------  --------------------------       ---------------- */
/* 1.0      2013-10-22  Added first version number       JMC              */
/* 1.0      2013-10-23  Moved Profile to diff file       JMC              */
/* 1.0-1    2013-10-28  Changed UNPRoTECT count in       JMC              */
/*                      prepReturnList as was                             */
/*                      causing seg faults                                */
/* 1.0-1    2013-10-28  Changed UNPRoTECT count in       JMC              */
/*                      prepReturnList as was                             */
/*                      causing seg faults                                */ 
/* 1.0-2    2013-10-29  Changed UNPRoTECT back in        JMC              */
/*                      prepReturnList, added                             */
/*                      UNPROTECT(1) in compare                           */ 
/* 1.0-3    2013-11-01  Calls to profile->compare        JMC              */
/*                      had rare and wildcard switces                     */
/*                      reversed.                                         */ 
/* 1.0-4    2013-11-04  Added some error checking       JMC               */
/*                      to prepReturnList                                 */
/* 1.0-5    2013-11-05  Changed compare                 JMC               */
/*                      arg lists and code                                */
/*                      to match so we don't guess that                   */
/*                      the input is correct                              */
/* 1.0-5    2013-11-05  Added trace code to             JMC               */
/*                      compare, mcompare to help                         */
/*                      debugging                                         */
/* 1.0-6    2018-03-22  Prototypes in header file.       MMA              */
/*                      Small adjustments.                                */


#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>

using namespace std;

//#undef length // remove due to Rcpp update?
 
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "compare-utils.h"

// [[Rcpp::export]]
Rcpp::List compare(const Rcpp::StringVector& DB, int numLoci, int bigHit, bool trace, int single,
                   bool useWildcard, bool useWildcardEffect, bool useRallele) {  
  if(trace){
    Rprintf("numLoci: %d\n", numLoci);
    Rprintf("bigHit: %d\n", bigHit);
    Rprintf("single: %d\n", single);
    Rprintf("useWildcard: %c\n", useWildcard ? 'T' : 'F');
    Rprintf("useWildcardEffect: %c\n", useWildcardEffect ? 'T' : 'F');
    Rprintf("useRallele: %c\n", useRallele ? 'T' : 'F');
  }
  
  vector<Profile*> vpProfiles;
  int nProfiles = DB.size();
  
  string strLine;
  string strID, strA1, strA2;
  
  int iProfiles = nProfiles;
  long unsigned comps = nProfiles*(nProfiles-1)/2;
  
  if (single > 0){ 
    iProfiles = single; // i only runs through 0...iProfiles-1
    comps = nProfiles * iProfiles; 
  }
  
  long unsigned stepper = 0;
  //long unsigned r = 0;
  
  Progress progress(comps, trace);
  

  // CONSTRUCT THE PROFILES VECTOR BY READING IN DATA FROM DB
  vpProfiles = readProfiles(DB, nProfiles, numLoci);
  
  unsigned long i, j;
  unsigned long m0, m1, m2, fm1, fm2;
  
  unsigned long nNumRows = useWildcardEffect ? 2 * numLoci + 1 : numLoci + 1;
  unsigned long m_size = nNumRows * nNumRows;
  
  Rcpp::IntegerVector m(m_size, 0);
//  PROTECT(m = allocVector(INTSXP, m_size));
/*  for(i = 0; i < m_size; i++){
    INTEGER(m)[i] = 0; 
  }*/
  
  vector<int> row1;
  vector<int> row2;
  vector<int> match;
  vector<int> partial;
  vector<int> fmatch;
  vector<int> fpartial;
  
  int ii;
  
  Profile *pProf1, *pProf2;
  
  // NEW!
  // the casts to unsigned long are mine (James) - I doubt it makes any difference
  for(i = 0; i < (unsigned long)iProfiles; i++){
    pProf1 = vpProfiles[i];
    if(single > 0){
      ii = single; 
    }
    else{
      ii = i+1;
    }
    
    
    for(j = ii; j < (unsigned long)nProfiles; j++){ // NEW! ends
    
    pProf2 = vpProfiles[j];
    
    m2 = 0;
    m1 = 0;
    fm2 = 0;
    fm1 = 0;
    
    pProf1->compare(pProf2, m2, m1, m0, fm2, fm1, useWildcardEffect, useRallele); // v1.0-3 Wildcard and Rare were reversed
    
    progress.increment();
    
    if(nProfiles >= 15 && trace){  // 15*14/2 = 105 > 100 whereas 14*13/2 = 91 < 100
      if (stepper > comps/100){ 
        if (Progress::check_abort()){
          Rcpp::stop("Aborted"); // FIXME: Return intermediate result?
        }
        stepper = 0;
      }
      stepper++;
    }
    
    //      m(m2,m1)++;
    if(useWildcardEffect){
      m[(m2 * 2 + m1) * (2 * numLoci + 1) + ( fm2 * 2 + fm1)]++;
    }
    else{
      m[(m2 + fm2) * (numLoci + 1)+( fm1 + m1)]++;
      /*
      if (m[2] > 7099) {      
        Rcpp::Rcout << "i = " << i << "; j = " << j << "; m[2] = " << m[2] << std::endl; 
        //Rcpp::print();
        //Rcpp::Rcout << m[2] << std::endl;
      }
      */
      
      if((m2 + fm2) >= (long unsigned)bigHit){
        // 	prof1.push_back(pProf1->m_strName);
        // 	prof2.push_back(pProf2->m_strName);
        row1.push_back(i + 1);
        row2.push_back(j + 1);
        match.push_back(m2);
        partial.push_back(m1);
        fmatch.push_back(fm2);
        fpartial.push_back(fm1);
      }
    }
    
    } // end for(j)
  } // end for(i)
  
  //  for(i=0;i<(numLoci+1)*(numLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);
  
  Rcpp::List rl = prepReturnList(m, row1, row2, match, partial, fmatch, fpartial);  
  
  return rl;
  
}

