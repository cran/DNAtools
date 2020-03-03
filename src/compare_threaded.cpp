#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>

using namespace std;

#include <Rcpp.h>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;
#include <tthread/fast_mutex.h>


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "compare-utils.h"




struct CompareWorker : public Worker {

  const size_t m_nProfiles;
  const int m_numLoci;
  
  const size_t m_bigHit;
  const int m_single;  
  const bool m_useWildcard;
  const bool m_useWildcardEffect;
  const bool m_useRallele;
  
  const vector<Profile*> &m_vpProfiles;

  unsigned long m_nNumRows;
  
  tthread::mutex m_mutex_out_m;
  tthread::mutex m_mutex_out_vectors;
  Rcpp::IntegerVector& out_m;
  
  vector<int>& out_row1;
  vector<int>& out_row2;
  vector<int>& out_match;
  vector<int>& out_partial;
  vector<int>& out_fmatch;
  vector<int>& out_fpartial;

  CompareWorker(const vector<Profile*> &vpProfiles, 
                size_t nProfiles,
                const int numLoci,

                const size_t bigHit,
                const int single,
                const bool useWildcard,
                const bool useWildcardEffect,
                const bool useRallele,
                
                Rcpp::IntegerVector &m,
                vector<int> &row1,
                vector<int> &row2,
                vector<int> &match,
                vector<int> &partial,
                vector<int> &fmatch,
                vector<int> &fpartial
                ) : 
                
    m_nProfiles(nProfiles), 
    m_numLoci(numLoci), 
    
    m_bigHit(bigHit), 
    m_single(single), 
    m_useWildcard(useWildcard), 
    m_useWildcardEffect(useWildcardEffect),
    m_useRallele(useRallele),
    
    m_vpProfiles(vpProfiles),
    
    m_nNumRows(m_useWildcardEffect ? 2 * m_numLoci + 1 : m_numLoci + 1),

    out_m(m),
    out_row1(row1),
    out_row2(row2),
    out_match(match),
    out_partial(partial),
    out_fmatch(fmatch),
    out_fpartial(fpartial) {
    
    if (out_m.size() != (m_nNumRows * m_nNumRows)) {
      Rcpp::stop("The m provided has the wrong size");
    }
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t ii, j;
    
    Profile *pProf1, *pProf2;
    
    unsigned long m0, m1, m2, fm1, fm2;
    // FIXME mikl: m0 ??
    
    //std::cout << "[" << begin << "; " << end << ")" << std::endl;
    
    // All have same size
    vector<int> local_row1;
    vector<int> local_row2;
    vector<int> local_match;
    vector<int> local_partial;
    vector<int> local_fmatch;
    vector<int> local_fpartial;
    
    vector<unsigned long> local_m_indices_increment;
    local_m_indices_increment.reserve(end - begin);

    // NEW!
    // the casts to unsigned long are mine (James) - I doubt it makes any difference
    for(std::size_t i = begin; i < end; i++){
      pProf1 = m_vpProfiles[i];
      if(m_single > 0){
        ii = m_single; 
      }
      else{
        ii = i+1;
      }


      for(j = ii; j < (unsigned long)m_nProfiles; j++){ // NEW! ends

      pProf2 = m_vpProfiles[j];

      m2 = 0;
      m1 = 0;
      fm2 = 0;
      fm1 = 0;

      pProf1->compare(pProf2, m2, m1, m0, fm2, fm1, m_useWildcardEffect, m_useRallele); // v1.0-3 Wildcard and Rare were reversed

      /*
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
      */
      
      //      m(m2,m1)++;
      if(m_useWildcardEffect){
        unsigned long idx = (m2 * 2 + m1) * (2 * m_numLoci + 1) + ( fm2 * 2 + fm1);
        local_m_indices_increment.push_back(idx);
      }
      else{
        unsigned long idx = (m2 + fm2) * (m_numLoci + 1)+( fm1 + m1);
        local_m_indices_increment.push_back(idx);
        
        /*
        if (m[2] > 7099) {      
          Rcpp::Rcout << "i = " << i << "; j = " << j << "; m[2] = " << m[2] << std::endl; 
          //Rcpp::print();
          //Rcpp::Rcout << m[2] << std::endl;
        }
        */
        
        if((m2 + fm2) >= (long unsigned)m_bigHit){
          local_row1.push_back(i + 1);
          local_row2.push_back(j + 1);
          local_match.push_back(m2);
          local_partial.push_back(m1);
          local_fmatch.push_back(fm2);
          local_fpartial.push_back(fm1);
        }
      }

      } // end for(j)
    } // end for(i)
    
    // Now move local containers' content to global result
    m_mutex_out_m.lock();
    for (auto idx : local_m_indices_increment) {
      out_m[idx]++;
    }
    m_mutex_out_m.unlock();    
    
    
    // All have same size
    m_mutex_out_vectors.lock();
    out_row1.reserve(out_row1.size() + local_row1.size());
    out_row2.reserve(out_row2.size() + local_row2.size());    
    out_match.reserve(out_match.size() + local_match.size());
    out_partial.reserve(out_partial.size() + local_partial.size());
    out_fmatch.reserve(out_fmatch.size() + local_fmatch.size());
    out_fpartial.reserve(out_fpartial.size() + local_fpartial.size());    
    for (size_t veci = 0; veci < local_row1.size(); ++veci) {
      out_row1.push_back( local_row1[veci] );
      out_row2.push_back( local_row2[veci] );
      out_match.push_back( local_match[veci] );
      out_partial.push_back( local_partial[veci] );
      out_fmatch.push_back( local_fmatch[veci] );
      out_fpartial.push_back( local_fpartial[veci] );    
    }
    m_mutex_out_vectors.unlock();
  }
};



// [[Rcpp::export]]
Rcpp::List compare_threaded(const Rcpp::StringVector& DB, int numLoci, int bigHit, bool trace, int single,
                            bool useWildcard, bool useWildcardEffect, bool useRallele) {  
                            
  if(trace){
    Rprintf("threaded\n");
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
  
  //long unsigned stepper = 0;
  
  Progress progress(comps, trace);

  // CONSTRUCT THE PROFILES VECTOR BY READING IN DATA FROM DB
  vpProfiles = readProfiles(DB, nProfiles, numLoci);
  
  unsigned long nNumRows = useWildcardEffect ? 2 * numLoci + 1 : numLoci + 1;
  unsigned long m_size = nNumRows * nNumRows;
  
  Rcpp::IntegerVector m(m_size, 0);
  
  vector<int> row1;
  vector<int> row2;
  vector<int> match;
  vector<int> partial;
  vector<int> fmatch;
  vector<int> fpartial;
  
  // NEW!
  // the casts to unsigned long are mine (James) - I doubt it makes any difference
  
  CompareWorker comp_work(vpProfiles, 
                          nProfiles, numLoci, 
                          bigHit, single,
                          useWildcard, useWildcard, useRallele,                                                    
                          m,
                          row1, row2, match, partial, fmatch, fpartial);
                            
  // FIXME: stack seem fragile to grain size
  //RcppParallel::parallelFor(0, (unsigned long)iProfiles, comp_work);
  RcppParallel::parallelFor(0, (unsigned long)iProfiles, comp_work, 1000);
  
  
  Rcpp::List rl = prepReturnList(m, row1, row2, match, partial, fmatch, fpartial);  
  
  return rl;
  
}

