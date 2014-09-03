/*########################################################################*/
/* File: mcompare.cpp                                                     */
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
/*                      UNPROTECT(1) in compare and                       */ 
/*                      UNPROTECT(1) in mcompare                          */ 
/* 1.0-3    2013-11-01  Calls to profile->compare        JMC              */
/*                      had rare and wildcard switces                     */
/*                      reversed.                                         */ 
/* 1.0-4    2013-11-04  Added some error checking       JMC               */
/*                      to prepReturnList                                 */
/* 1.0-5    2013-11-05  Uncommented UNPROTECT in        JMC               */
/*                      mcompare                                          */
/* 1.0-5    2013-11-05  Changed compare, mcompare       JMC               */
/*                      arg lists and code                                */
/*                      to match so we don't guess that                   */
/*                      the input is correct                              */
/* 1.0-5    2013-11-05  Added trace code to             JMC               */
/*                      compare, mcompare to help                         */
/*                      debugging                                         */




#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rcpp.h>

// JMC: header for interrupts
#include <R_ext/Utils.h>

// JMC: removed Profile class declaration and implementation to separate files for maintainability
#include "profile.h"

// Mikkel:
#ifdef MACOS
#define DBCOMPARE_MULTICORE
#include <sys/param.h>
#include <sys/sysctl.h>
#include <unistd.h>
#elif __linux__
#define DBCOMPARE_MULTICORE
#include <unistd.h>
#else // Windows, Solaris, BSD etc.
// Nothing here
#endif

// From http://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine
int getNumCores() {
  #ifdef MACOS
  int nm[2];
  size_t len = 4;
  uint32_t count;
  
  nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
  sysctl(nm, 2, &count, &len, NULL, 0);
  
  if(count < 1) {
    nm[1] = HW_NCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);
    if(count < 1) { count = 1; }
  }
  
  return count;
  #elif __linux__
  return sysconf(_SC_NPROCESSORS_ONLN);
  #else
  return 1;
  #endif
}

void my_usleep() {
  #if defined MACOS || defined __linux__ 
  usleep(500000);
  #endif
} 

vector<Profile *> readProfiles(Rcpp::StringVector &DB, int nProfiles, int nLoci){
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

#ifdef DBCOMPARE_MULTICORE
class thread_data {
  public:
  pthread_mutex_t *change_lock;
  pthread_mutex_t *i_queue_lock;
  int *i;
  int id;
  int nProfiles;
  int Profiles;
  int nLoci;
  int hit;
  bool bTrace;
  int nSingle;
  bool bWildcard; // wildcard modification
  bool bWildcardEffect; // wildcard modification
  bool bRallele;
  long unsigned* stepper;
  vector<Profile*> *vpProfiles;
  Rcpp::IntegerVector m;
  Rcpp::IntegerVector M;
  vector<int> *row1;
  vector<int> *row2;
  vector<int> *match;
  vector<int> *partial;
  vector<int> *fmatch; // wildcard modification
  vector<int> *fpartial; // wildcard modification
};
#endif

#ifdef DBCOMPARE_MULTICORE
void *compare_profiles(void *threadarg) {
  thread_data *my_data;
  my_data = (thread_data*) threadarg;
  
  pthread_mutex_t* change_lock = my_data->change_lock;
  pthread_mutex_t* i_queue_lock = my_data->i_queue_lock;
  
  //HER:	int id = my_data->id;
  int nProfiles = my_data->nProfiles;
  int nLoci = my_data->nLoci;
  int hit = my_data->hit;
  //HER:	int trace = my_data->trace;
  int nSingle = my_data->nSingle;
  bool bWildcard = my_data->bWildcard;  // wildcard modification
  bool bWildcardEffect = my_data->bWildcardEffect;  // wildcard modification
  bool bRallele = my_data->bRallele;
  long unsigned *stepper = my_data->stepper;
  vector<Profile*> *vpProfiles = my_data->vpProfiles;
  Rcpp::IntegerVector m = my_data->m;
  vector<int> *row1 = my_data->row1;
  vector<int> *row2 = my_data->row2;
  vector<int> *match = my_data->match;
  vector<int> *partial = my_data->partial;
  vector<int> *fmatch = my_data->fmatch; // wildcard modification
  vector<int> *fpartial = my_data->fpartial; // wildcard modification
  
  vector<int> local_row1;
  vector<int> local_row2;
  vector<int> local_match;
  vector<int> local_partial;
  vector<int> local_fmatch; // wildcard modification
  vector<int> local_fpartial; // wildcard modification
  
  long unsigned i, j, k, m0, m1, m2, fm1, fm2; // wildcard modification
  int nA1, nA2, nB1, nB2;
  
  int nNumRows = bWildcardEffect ? 2 * nLoci + 1 : nLoci + 1;
  int m_size = nNumRows * nNumRows; // wildcard modification
  
  Rcpp::IntegerVector local_m(m_size, 0);
  
  Profile *pProf1, *pProf2;
  
  // Mikkel:
  
  for (;;) { // we break when there's no more i's left
  
  //// SOMETHING HERE ABOUT single>0 AND HOW TO SET i
  pthread_mutex_lock(i_queue_lock);
  i = *my_data->i;
  *my_data->i += 1;
  pthread_mutex_unlock(i_queue_lock);
  
  if (i >= (long unsigned)nProfiles) {
    break;
  }
  
  //Rprintf("Thread %i calculates i = %i\n", id, i);
  
  pProf1 = (*vpProfiles)[i];
  
  for(j = i + 1; j < (long unsigned)nProfiles; j++){
    pProf2 = (*vpProfiles)[j];
    
    m2 = 0;
    m1 = 0;
    m0 = 0;
    fm2 = 0; // wildcard modification
    fm1 = 0; // wildcard modification
    
    pProf1->compare(pProf2, m2, m1, m0, fm2, fm1, bWildcardEffect, bRallele); // 1.0-3 Wildcard and rare were reversed
    
    // thread_m(m2,m1)++;
    if(bWildcardEffect){
      local_m[(m2 * 2 + m1) * (2 * nLoci + 1) + ( fm2 * 2 + fm1)] ++;
    }
    else{
      local_m[(m2 + fm2) * (nLoci + 1) + (m1 + fm1)] ++;
      
      if((m2 + fm2) >= (long unsigned)hit){
        //Rprintf("hit = %i, m2 = %i\n", hit, m2);
        //	 prof1.push_back(pProf1->m_strName);
        //	 prof2.push_back(pProf2->m_strName);
        local_row1.push_back(i+1);
        local_row2.push_back(j+1);
        local_match.push_back(m2);
        local_partial.push_back(m1);
        local_fmatch.push_back(fm2);
        local_fpartial.push_back(fm1);
      }
    }
    
  } // end for(j)
  
  pthread_mutex_lock(change_lock); // lock because we need to write to the shared data
  *stepper += nProfiles - (i+1);
  pthread_mutex_unlock(change_lock); // unlock again
  }
  
  if(!bWildcardEffect){
    int n1 = local_row1.size();
    
    pthread_mutex_lock(change_lock); // lock because we need to write to the shared data
    for(j = 0; j < (long unsigned)n1; j++) {
      row1->push_back(local_row1[j]);
      row2->push_back(local_row2[j]);
      match->push_back(local_match[j]);
      partial->push_back(local_partial[j]);
      fmatch->push_back(local_fmatch[j]);
      fpartial->push_back(local_fpartial[j]);
    }
  }
  
  for(j = 0; j < (long unsigned)m_size; j++){ 
    // m(m2,m1) + thread_m(m2,m1)
    m[j] += local_m[j];
  }
  pthread_mutex_unlock(change_lock); // unlock again
   
  pthread_exit(NULL);
}
#endif

//// SINGLE CORE CALL
// [[Rcpp::export]]
Rcpp::List compare(SEXP db, SEXP numLoci, SEXP bigHit, SEXP trace, SEXP single,
                   SEXP useWildcard, SEXP useWildcardEffect, SEXP useRallele) {
  
  Rcpp::StringVector DB(db);
  
  int nLoci = Rcpp::as<int>(numLoci);
  int hit = Rcpp::as<int>(bigHit);
  bool bTrace = (bool)Rcpp::as<int>(trace);
  int nSingle = Rcpp::as<int>(single);
  bool bWildcard = (bool)Rcpp::as<int>(useWildcard);
  bool bWildcardEffect = (bool)Rcpp::as<int>(useWildcardEffect);
  bool bRallele = (bool)Rcpp::as<int>(useRallele);
  
  if(bTrace){
    Rprintf("Single threaded call\n");
    Rprintf("nLoci: %d\n", nLoci);
    Rprintf("hit: %d\n", hit);
    Rprintf("nSingle: %d\n", nSingle);
    Rprintf("bWildcard: %c\n", bWildcard ? 'T' : 'F');
    Rprintf("bWildcardEffect: %c\n", bWildcardEffect ? 'T' : 'F');
    Rprintf("bRallele: %c\n", bRallele ? 'T' : 'F');
  }
  
  
  vector<Profile*> vpProfiles;
  int nProfiles = DB.size();
  
  string strLine;
  string strID, strA1, strA2;
  
  int iProfiles = nProfiles;
  long unsigned comps = nProfiles*(nProfiles-1)/2;
  
  if(nSingle > 0){ 
    iProfiles = nSingle; // i only runs through 0...iProfiles-1
    comps = nProfiles * iProfiles; 
  }
  
  long unsigned stepper = 0;
  long unsigned r = 0;
  
  if(bTrace){
    Rprintf("Progress:\n");
    for(r=0; r < 101; r = r + 5){
      if((r % 10) == 0){ 
        Rprintf("%d%%",r); 
      }else{
        Rprintf("       ");
      }
    }
    Rprintf("\n");
    
    for(r = 0; r < 101; r++){
      if((r % 10) == 0){
        Rprintf(":");
      }else{
        Rprintf(".");
      }
    }
    Rprintf("\n|");
  }
  
  // CONSTRUCT THE PROFILES VECTOR BY READING IN DATA FROM DB
  vpProfiles = readProfiles(DB, nProfiles, nLoci);
  
  unsigned long i, j;
  unsigned long m0, m1, m2, fm1, fm2;
  
  unsigned long nNumRows = bWildcardEffect ? 2 * nLoci + 1 : nLoci + 1;
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
    if(nSingle > 0){
      ii = nSingle; 
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
    
    pProf1->compare(pProf2, m2, m1, m0, fm2, fm1, bWildcardEffect, bRallele); // v1.0-3 Wildcard and Rare were reversed
    
    
    if(nProfiles >= 15 && bTrace){  // 15*14/2 = 105 > 100 whereas 14*13/2 = 91 < 100
      if(stepper > comps/100){ 
        Rprintf("=");
        R_CheckUserInterrupt();
        stepper = 0;
      }
      stepper++;
    }
    
    //      m(m2,m1)++;
    if(bWildcardEffect){
      m[(m2 * 2 + m1) * (2 * nLoci + 1) + ( fm2 * 2 + fm1)]++;
    }
    else{
      m[(m2 + fm2) * (nLoci + 1)+( fm1 + m1)]++;
      
      if((m2 + fm2) >= (long unsigned)hit){
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
  
  if(bTrace){
    if(nProfiles < 15){ 
      for(r = 0; r < 99; r++)
      Rprintf("=");
    }
    Rprintf("|\n");
  }
  
  //  for(i=0;i<(nLoci+1)*(nLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);
  
  Rcpp::List rl = prepReturnList(m, row1, row2, match, partial, fmatch, fpartial);  
  
  return rl;
  
}

// [[Rcpp::export]]
Rcpp::IntegerVector score(SEXP prof1, SEXP prof2, SEXP numLoci, SEXP useWildCard, SEXP useRareAllele){
  
  int *pnProf1 = &(Rcpp::as<std::vector<int> >(prof1)[0]); // this should replace INTEGER(prof1)
  int *pnProf2 = &(Rcpp::as<std::vector<int> >(prof2)[0]);
  
  int nLoci = Rcpp::as<int>(numLoci);
  bool bWildCard = (bool)(Rcpp::as<int>(useWildCard));
  bool bRareAllele = (bool)(Rcpp::as<int>(useRareAllele));
  
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

//// MULTI CORE CALL
// [[Rcpp::export]]
Rcpp::List mcompare(SEXP db, SEXP numLoci, SEXP bigHit, SEXP trace, SEXP single, SEXP numThreads, 
              SEXP useWildcard, SEXP useWildcardEffect, SEXP useRallele) { 
  #ifndef DBCOMPARE_MULTICORE
  return compare(db, numLoci, bigHit, trace, single, useWildcard, useWildcardEffect, useRallele);
  #else
  
  Rcpp::StringVector DB(db);
  
  int nLoci = Rcpp::as<int>(numLoci);
  int hit = Rcpp::as<int>(bigHit);
  bool bTrace = (bool)Rcpp::as<int>(trace);
  int nSingle = Rcpp::as<int>(single);
  int threads = Rcpp::as<int>(numThreads);
  bool bWildcard = (bool)Rcpp::as<int>(useWildcard);
  bool bWildcardEffect = (bool)Rcpp::as<int>(useWildcardEffect);
  bool bRallele = (bool)Rcpp::as<int>(useRallele);
  
  if(bTrace){
    Rprintf("Multithreaded call\n");
    Rprintf("nLoci: %d\n", nLoci);
    Rprintf("hit: %d\n", hit);
    Rprintf("nSingle: %d\n", nSingle);
    Rprintf("threads: %d\n", threads);
    Rprintf("bWildcard: %c\n", bWildcard ? 'T' : 'F');
    Rprintf("bWildcardEffect: %c\n", bWildcardEffect ? 'T' : 'F');
    Rprintf("bRallele: %c\n", bRallele ? 'T' : 'F');
  }
  
  if (threads == 1) {
    return compare(db, numLoci, bigHit, trace, single, useWildcard, useWildcardEffect, useRallele);
  }
  
  int number_of_threads = getNumCores();
  //number_of_threads = 1;
  
  if (threads > 0 && threads <= number_of_threads) { // fewer threads than cores requested
  number_of_threads = threads;
  }
  
  vector<Profile*> vpProfiles;
  int nProfiles = DB.size();
  
  long unsigned comps = nProfiles * (nProfiles - 1) / 2;
  
  long unsigned stepper = 0;
  int r;
  
  if(bTrace){
    Rprintf("Using %i threads\n", number_of_threads);
    
    Rprintf("Progress:\n");
    for(r = 0; r < 101; r = r+5){
      if((r % 10) == 0){ Rprintf("%d%%",r); }
      else{ Rprintf("       "); }
    }
    Rprintf("\n");
    
    for(r=0;r<101;r++){
      if((r%10) == 0){ Rprintf(":"); }
      else{ Rprintf("."); }
    }
    Rprintf("\n|");
  }
  
  // CONSTRUCT THE PROFILES VECTOR BY READING IN DATA FROM DB
  vpProfiles = readProfiles(DB, nProfiles, nLoci);
  
  int i, j;
  
  int nNumRows = bWildcardEffect ? 2 * nLoci + 1 : nLoci + 1;
  int m_size = nNumRows * nNumRows;
  Rcpp::IntegerVector m(m_size, 0);
  
  vector<int> row1;
  vector<int> row2;
  vector<int> match;
  vector<int> partial;
  vector<int> fmatch;
  vector<int> fpartial;
  
  vector<thread_data> thread_data_array;
  
  int rc, status;
  pthread_mutex_t change_lock;
  pthread_mutex_t i_queue_lock;
  
  if (pthread_mutex_init(&change_lock, NULL)) {
    // REprintf("Could not initialize change_lock mutex, aborting.\n");
    // fprintf(stderr, "Could not initialize change_lock mutex, aborting.\n");
    error("Could not initialize change_lock mutex, aborting.\n"); 
    // exit(1);
  }
  
  if (pthread_mutex_init(&i_queue_lock, NULL)) {
    // REprintf("Could not initialize i_queue_lock mutex, aborting.\n");
    // fprintf(stderr, "Could not initialize i_queue_lock mutex, aborting.\n");
    error("Could not initialize i_queue_lock mutex, aborting.\n"); 
    // exit(1);
  }
  
  vector<pthread_t> threads_container;
  thread_data *mthread_data;
  
  int i_row = 0;
  
  for(j = 0; j < number_of_threads; j++){
    //Rprintf("Thread %i initiating... ", j);
    mthread_data = new thread_data();
    mthread_data->change_lock = &change_lock;
    mthread_data->i_queue_lock = &i_queue_lock;
    mthread_data->i = &i_row;
    mthread_data->id = j;
    mthread_data->nProfiles = nProfiles;
    mthread_data->nLoci = nLoci;
    mthread_data->hit = hit;
    mthread_data->bTrace = bTrace;
    mthread_data->nSingle = nSingle;
    mthread_data->bWildcard = bWildcard;
    mthread_data->bWildcardEffect = bWildcardEffect;
    mthread_data->bRallele = bRallele;
    mthread_data->stepper = &stepper;
    mthread_data->vpProfiles = &vpProfiles;
    mthread_data->m = m;
    mthread_data->row1 = &row1;
    mthread_data->row2 = &row2;
    mthread_data->match = &match;
    mthread_data->partial = &partial;
    mthread_data->fmatch = &fmatch;
    mthread_data->fpartial = &fpartial;
    thread_data_array.push_back(*mthread_data);
    
    pthread_t thr;
    rc = pthread_create(&thr, NULL, compare_profiles, (void *) mthread_data);
    threads_container.push_back(thr);
    
    if (rc) {
      // REprintf("Error in thread creation: %i.\n", rc);
      // fprintf(stderr, "Error in thread creation: %i.\n", rc);
      error("Error in thread creation");
      // exit(-1);
    }
  }
  
  int progress_set = 0;
  long unsigned comp_div = comps/100;
  
  if(nProfiles >= 15 && bTrace){
    for(;;) {
      if(progress_set >= 99) {
        break;
      }
      
      if (i_row >= nProfiles) { // no need to lock when only reading i
      break;
      } 
      
      // 15*14/2 = 105 > 100 whereas 14*13/2 = 91 < 100
      if(stepper > comp_div){			 
        while (stepper > comp_div){
          Rprintf("=");
          progress_set++;
          
          pthread_mutex_lock(&change_lock);
          stepper -= comp_div;
          pthread_mutex_unlock(&change_lock);
        }
      }
      
      my_usleep();
    }
  }
  
  // print the missing progress bar, the above is only approx due to int div
  
  if(bTrace){
    for(j = progress_set; j < 99; j++){
      Rprintf("=");
    }
  }
  
  for(j = 0; j < number_of_threads; j++){
    rc = pthread_join(threads_container[j], (void **)&status);
    
    if (rc) {
      // REprintf("Error in thread joining: %i.\n", rc); 
      // fprintf(stderr, "Error in thread joining: %i.\n", rc);
      error("Error in thread joining"); 
      // exit(-1);
    }
  }
  
  if(bTrace){
    if(nProfiles < 15){ 
      for(r = 0; r < 99; r++) 
      Rprintf("=");
    }
    Rprintf("|\n");
  }
  
  //	for(i=0;i<(nLoci+1)*(nLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);
  Rcpp::List rl = prepReturnList(m, row1, row2, match, partial, fmatch, fpartial);  
  
  return rl;
  
  #endif
}

