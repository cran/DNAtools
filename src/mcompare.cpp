
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
// Mikkel:
#include <queue>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

// Mikkel:
#ifdef MACOS
#define DBCOMPARE_MULTICORE
#include <sys/param.h>
#include <sys/sysctl.h>
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

extern "C" { // begin extern
  
  // From Rcpp
  class StringVector {
  public:
    StringVector(SEXP vec);
    ~StringVector() {
      delete [] v;
    }
    inline std::string& operator()(int i) {
      return v[i];
    }
	int size() { return length; }
  private:
    std::string *v;
    int length;
  };
  
  StringVector::StringVector(SEXP vec) {
    int i;
    int len = length(vec);
    v = new std::string[len];
    for(i = 0; i < len; i++)
      v[i] = std::string(CHAR(STRING_ELT(vec,i)));
    length = len;
  }
  
  using namespace std;
  
  class Profile{
  public:
    int **m_ppdProfile;
    int m_nLoci;
    
  public:
    Profile(void){
      m_ppdProfile = NULL;
    };
    
    Profile(int nLoci){
      m_nLoci = nLoci;
      m_ppdProfile = new int*[nLoci];
      int nLoc;
      for(nLoc=0;nLoc<nLoci;nLoc++)
	m_ppdProfile[nLoc] = new int[2];
    };
    
    Profile(const Profile &p){
      m_nLoci = p.m_nLoci;
      m_ppdProfile = new int*[m_nLoci];
      int nLoc;
      for(nLoc=0;nLoc<m_nLoci;nLoc++){
	m_ppdProfile[nLoc] = new int[2];
	m_ppdProfile[nLoc][0] = p.m_ppdProfile[nLoc][0];
	m_ppdProfile[nLoc][1] = p.m_ppdProfile[nLoc][1];
      }
    };
    
    const Profile& operator=(const Profile &p){
      m_nLoci = p.m_nLoci;
      m_ppdProfile = new int*[m_nLoci];
      int nLoc;
      for(nLoc=0;nLoc<m_nLoci;nLoc++){
	m_ppdProfile[nLoc] = new int[2];
	m_ppdProfile[nLoc][0] = p.m_ppdProfile[nLoc][0];
	m_ppdProfile[nLoc][1] = p.m_ppdProfile[nLoc][1];
      }
      
      return *this;
    };
    
    vector<int> Compare(const Profile &p){
      vector<int> vResult;
      return vResult;
    };
  };

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
    int trace;
    int single;
    int wildcard; // wildcard modification
    int wildcardEffect; // wildcard modification
    long unsigned* stepper;
    vector<Profile*> *vpProfiles;
    SEXP m;
    SEXP M;
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
    int single = my_data->single;
    int wildcard = my_data->wildcard;  // wildcard modification
    int wildcardEffect = my_data->wildcardEffect;  // wildcard modification
    long unsigned *stepper = my_data->stepper;
    vector<Profile*> *vpProfiles = my_data->vpProfiles;
    SEXP m = my_data->m;
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
    
    int i, j, k, m1, m2, fm1, fm2; // wildcard modification
    int dA1, dA2, dB1, dB2;
    
    int m_size = ((1+wildcardEffect)*nLoci+1)*((1+wildcardEffect)*nLoci+1); // wildcard modification
    
    SEXP local_m;
    PROTECT(local_m = allocVector(INTSXP, m_size));
    
    for(j=0;j<m_size;j++){ INTEGER(local_m)[j] = 0; }
    
    Profile *pProf1, *pProf2;
    
    // Mikkel:
    
    for (;;) { // we break when there's no more i's left
      
      //// SOMETHING HERE ABOUT single>0 AND HOW TO SET i
      pthread_mutex_lock(i_queue_lock);
      i = *my_data->i;
      *my_data->i += 1;
      pthread_mutex_unlock(i_queue_lock);
      
      if (i >= nProfiles) {
	break;
      }
      
      //Rprintf("Thread %i calculates i = %i\n", id, i);
      
      pProf1 = (*vpProfiles)[i];
      
      for(j=i+1;j<nProfiles;j++){
	pProf2 = (*vpProfiles)[j];

	m2 = 0;
	m1 = 0;
	fm2 = 0; // wildcard modification
	fm1 = 0; // wildcard modification
	
	for(k=0;k<nLoci;k++){
	  dA1 = pProf1->m_ppdProfile[k][0];
	  dA2 = pProf1->m_ppdProfile[k][1];
	  dB1 = pProf2->m_ppdProfile[k][0];
	  dB2 = pProf2->m_ppdProfile[k][1];
	  
	  // if(wildcard==1){ // Overwrites profile if hom (aa) to get aF
	  //   if(dA1==dA2) dA1 = 0;
	  //   if(dB1==dB2) dB1 = 0;
	  // }
	  
	  if(dA1==dA2){ // Profile 1 is hom
	    if(dA1==0){ // Both A alleles are wildcards "F"
	      fm2++; // wildcard modification
	    }
	    else{ // No wildcards in profile 1
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++; // Both B alleles are wildcards "F"
		else if(dA1==dB1){
		  m2++; // A genuine match
		  if(wildcardEffect==1) fm2++; // would be aF,aF  // wildcard modification
		}
		else if(wildcardEffect==1) fm1++; // if analysing effect of F  // wildcard modification
	      }
	      else{ // Profile 2 is het
		if(dB1==0){ // Profile 2 is bF
		  if(dA1==dB2) fm2++;
		  else fm1++;
		}
		else{ // Profile 2 has no wildcards
		  if(dA1==dB1 || dA1==dB2){
		    m1++;
		    if(wildcardEffect==1) fm2++; // would be aF,ab
		  }
		  else if(wildcardEffect==1) fm1++;
		}
	      }
	    }
	  }
	  else{ // Profile 1 is het
	    if(dA1==0){ // Profile 1 is aF
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++;
		else if(dA2==dB1) fm2++;
		else fm1++;
	      }
	      else{ // Profile 2 is het
		if(dB1==0){
		  if(dA2==dB2) fm2++;
		  else fm1++;
		}
		else{ // Profile 2 has no wildcards
		  if(dA2==dB1 || dA2==dB2) fm2++;
		  else fm1++;
		}
	      }
	    }
	    else{ // profile 1 has no wildcards
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++;
		else if(dA1==dB1 || dA2==dB1){ 
		  m1++;
		  if(wildcardEffect==1) fm2++;
		}
		else if(wildcardEffect==1) fm1++;
	      }
	      else{ // Profile 2 is het
		if(dB1==0){
		  if(dA1==dB2 || dA2==dB2) fm2++;
		  else fm1++;
		}
		else{ // None of the profiles has wildcards
		  if(dA1==dB1 && dA2==dB2){
		    m2++;
		    if(wildcardEffect==1) fm2++;
		  }
		  else if((dA1==dB1 && dA2!=dB2) || (dA1==dB2 && dA2!=dB1) || (dA1!=dB1 && dA2==dB2) || (dA1!=dB2 && dA2==dB1)){
		    m1++;
		    if(wildcardEffect==1) fm1++;
		  }
		}
	      }
	    }
	  }
	  
	}
	
	// thread_m(m2,m1)++;
	if(wildcardEffect==1){
	  INTEGER(local_m)[(m2*2+m1)*(2*nLoci+1) + (fm2*2+fm1)] ++;
	}
	else{
	  INTEGER(local_m)[(m2+fm2)*(nLoci+1) + (m1+fm1)] ++;
	  
	  if((m2+fm2)>=hit){
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
    
    if(wildcardEffect==0){
      int n1 = local_row1.size();
      
      pthread_mutex_lock(change_lock); // lock because we need to write to the shared data
      for(j=0;j<n1;j++) {
	row1->push_back(local_row1[j]);
	row2->push_back(local_row2[j]);
	match->push_back(local_match[j]);
	partial->push_back(local_partial[j]);
	fmatch->push_back(local_fmatch[j]);
	fpartial->push_back(local_fpartial[j]);
      }
    }
    
    for(j=0;j<m_size;j++){ 
      // m(m2,m1) + thread_m(m2,m1)
      INTEGER(m)[j] += INTEGER(local_m)[j];
    }
    pthread_mutex_unlock(change_lock); // unlock again
    
    UNPROTECT(1); // local_m
    
    pthread_exit(NULL);
  }
#endif
  
  //// SINGLE CORE CALL
  
  SEXP compare(SEXP db, SEXP param) {
    
    StringVector DB(db);
    
    int nLoci = INTEGER(param)[0];
    int hit = INTEGER(param)[1];
    int trace = INTEGER(param)[2];
    int single = INTEGER(param)[3];
    int wildcard = INTEGER(param)[4];
    int wildcardEffect = INTEGER(param)[5];
    
    vector<Profile*> vpProfiles;
    int nProfiles = DB.size();
    
    string strLine;
    
    int nLoc;
    string strID, strA1, strA2;
    size_t nPos;
    Profile *pProfile;
    char cDelim = '\t';
    
    int iProfiles = nProfiles;
    long unsigned comps = nProfiles*(nProfiles-1)/2;
    
    if(single>0){ 
      iProfiles = single; // i only runs through 0...iProfiles-1
      comps = nProfiles*iProfiles; 
    }
    
    long unsigned stepper = 0;
    int r;
    
    if(trace==1){
      Rprintf("Progress:\n");
      for(r=0;r<101;r = r+5){
	if((r%10) == 0){ Rprintf("%d%%",r); }
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
    int t=0;
    while(t<nProfiles){
      strLine = DB(t);
      pProfile = new Profile(nLoci);
      for(nLoc=0;nLoc<nLoci;nLoc++){
	nPos = strLine.find(cDelim);
	strA1 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
	nPos = strLine.find(cDelim);
	strA2 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
	pProfile->m_ppdProfile[nLoc][0] = atoi(strA1.c_str());
	pProfile->m_ppdProfile[nLoc][1] = atoi(strA2.c_str());
      }
      vpProfiles.push_back(pProfile);
      t++;
    }
    
    int i, j, k, m1, m2, fm1, fm2;
    
    int m_size = ((1+wildcardEffect)*nLoci+1)*((1+wildcardEffect)*nLoci+1);
    
    SEXP m;
    PROTECT(m = allocVector(INTSXP, m_size));
    for(i=0;i<m_size;i++){ INTEGER(m)[i] = 0; }
    
    vector<int> row1;
    vector<int> row2;
    vector<int> match;
    vector<int> partial;
    vector<int> fmatch;
    vector<int> fpartial;
    
    int dA1, dA2, dB1, dB2, ii;
    
    Profile *pProf1, *pProf2;
    
    // NEW!
    for(i=0;i<iProfiles;i++){
      pProf1 = vpProfiles[i];
      if(single>0){ ii = single; }
      else{ ii = i+1; }
      for(j=ii;j<nProfiles;j++){ // NEW! ends
	
	pProf2 = vpProfiles[j];
	
	m2 = 0;
	m1 = 0;
	fm2 = 0;
	fm1 = 0;
	
	for(k=0;k<nLoci;k++){
	  dA1 = pProf1->m_ppdProfile[k][0];
	  dA2 = pProf1->m_ppdProfile[k][1];
	  dB1 = pProf2->m_ppdProfile[k][0];
	  dB2 = pProf2->m_ppdProfile[k][1];
	  
	  // if(wildcard==1){ // Overwrites the profiles - if homs (aa) we get aF
	  //   if(dA1==dA2) dA1 = 0;
	  //   if(dB1==dB2) dB1 = 0;
	  // }
	  
	  if(dA1==dA2){ // Profile 1 is hom
	    if(dA1==0){ // Both A alleles are wildcards "F"
	      fm2++;
	    }
	    else{ // No wildcards in profile 1
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++; // Both B alleles are wildcards "F"
		else if(dA1==dB1){
		  m2++; // A genuine match
		  if(wildcardEffect==1) fm2++; // would be aF,aF
		}
		else if(wildcardEffect==1) fm1++; // if analysing effect of F 
	      }
	      else{ // Profile 2 is het
		if(dB1==0){ // Profile 2 is bF
		  if(dA1==dB2) fm2++;
		  else fm1++;
		}
		else{ // Profile 2 has no wildcards
		  if(dA1==dB1 || dA1==dB2){
		    m1++;
		    if(wildcardEffect==1) fm2++; // would be aF,ab
		  }
		  else if(wildcardEffect==1) fm1++;
		}
	      }
	    }
	  }
	  else{ // Profile 1 is het
	    if(dA1==0){ // Profile 1 is aF
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++;
		else if(dA2==dB1) fm2++;
		else fm1++;
	      }
	      else{ // Profile 2 is het
		if(dB1==0){
		  if(dA2==dB2) fm2++;
		  else fm1++;
		}
		else{ // Profile 2 has no wildcards
		  if(dA2==dB1 || dA2==dB2) fm2++;
		  else fm1++;
		}
	      }
	    }
	    else{ // profile 1 has no wildcards
	      if(dB1==dB2){ // Profile 2 is hom
		if(dB1==0) fm2++;
		else if(dA1==dB1 || dA2==dB1){ 
		  m1++;
		  if(wildcardEffect==1) fm2++;
		}
		else if(wildcardEffect==1) fm1++;
	      }
	      else{ // Profile 2 is het
		if(dB1==0){
		  if(dA1==dB2 || dA2==dB2) fm2++;
		  else fm1++;
		}
		else{ // None of the profiles has wildcards
		  if(dA1==dB1 && dA2==dB2){
		    m2++;
		    if(wildcardEffect==1) fm2++;
		  }
		  else if((dA1==dB1 && dA2!=dB2) || (dA1==dB2 && dA2!=dB1) || (dA1!=dB1 && dA2==dB2) || (dA1!=dB2 && dA2==dB1)){
		    m1++;
		    if(wildcardEffect==1) fm1++;
		  }
		}
	      }
	    }
	  }
	  
	}
	
	if(nProfiles>=15 && trace==1){  // 15*14/2 = 105 > 100 whereas 14*13/2 = 91 < 100
	  if(stepper > comps/100){ 
	    Rprintf("=");
	    stepper = 0;
	  }
	  stepper++;
	}
	
	//      m(m2,m1)++;
	if(wildcardEffect==1){
	  INTEGER(m)[(m2*2+m1)*(2*nLoci+1) + (fm2*2+fm1)] ++;
	}
	else{
	  INTEGER(m)[(m2+fm2)*(nLoci+1)+(fm1+m1)] ++;
	  
	  if((m2+fm2)>=hit){
	    // 	prof1.push_back(pProf1->m_strName);
	    // 	prof2.push_back(pProf2->m_strName);
	    row1.push_back(i+1);
	    row2.push_back(j+1);
	    match.push_back(m2);
	    partial.push_back(m1);
	    fmatch.push_back(fm2);
	    fpartial.push_back(fm1);
	  }
	}
	
      } // end for(j)
    } // end for(i)
    
    if(trace==1){
      if(nProfiles<15){ 
	for(r=0;r<99;r++) Rprintf("=");
      }
      Rprintf("|\n");
    }
    
    //  for(i=0;i<(nLoci+1)*(nLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);
    
    SEXP r1, r2, mm, pp, fmm, fpp, list_names, list;
    PROTECT(r1 = allocVector(INTSXP, row1.size()));
    PROTECT(r2 = allocVector(INTSXP, row2.size()));
    PROTECT(mm = allocVector(INTSXP, match.size()));
    PROTECT(pp = allocVector(INTSXP, partial.size()));
    PROTECT(fmm = allocVector(INTSXP, fmatch.size()));
    PROTECT(fpp = allocVector(INTSXP, fpartial.size()));
    int matchlength = row1.size();
    for(i=0;i<matchlength;i++){
      INTEGER(r1)[i] = row1[i]; 
      INTEGER(r2)[i] = row2[i]; 
      INTEGER(mm)[i] = match[i]; 
      INTEGER(pp)[i] = partial[i]; 
      INTEGER(fmm)[i] = fmatch[i]; 
      INTEGER(fpp)[i] = fpartial[i]; 
    }
    
    PROTECT(list_names = allocVector(STRSXP,7));
    const char *names[7] = {"M","row1","row2","matches","partial","fmatches","fpartial"};
    for(i=0;i<7;i++) SET_STRING_ELT(list_names, i ,mkChar(names[i]));
    
    PROTECT(list = allocVector(VECSXP,7));
    SET_VECTOR_ELT(list, 0, m);
    SET_VECTOR_ELT(list, 1, r1);
    SET_VECTOR_ELT(list, 2, r2);
    SET_VECTOR_ELT(list, 3, mm);
    SET_VECTOR_ELT(list, 4, pp);
    SET_VECTOR_ELT(list, 5, fmm);
    SET_VECTOR_ELT(list, 6, fpp);
    setAttrib(list,R_NamesSymbol, list_names);
    
    UNPROTECT(9);
    
    return list;
    
  }
  
  //// MULTI CORE CALL
  
  SEXP mcompare(SEXP db, SEXP param) { 
#ifndef DBCOMPARE_MULTICORE
    return compare(db, param);
#else
    
    StringVector DB(db);
    
    int nLoci = INTEGER(param)[0];
    int hit = INTEGER(param)[1];
    int trace = INTEGER(param)[2];
    int single = INTEGER(param)[3];
    int threads = INTEGER(param)[4];
    int wildcard = INTEGER(param)[5];
    int wildcardEffect = INTEGER(param)[6];
    
    if (threads == 1) {
      return compare(db, param);
    }
    
    int number_of_threads = getNumCores();
    //number_of_threads = 1;
    
    if (threads > 0 && threads <= number_of_threads) { // fewer threads than cores requested
      number_of_threads = threads;
    }
    
    vector<Profile*> vpProfiles;
    int nProfiles = DB.size();
    
    string strLine;
    
    int nLoc;
    string strID, strA1, strA2;
    size_t nPos;
    Profile *pProfile;
    char cDelim = '\t';
    
    long unsigned comps = nProfiles*(nProfiles-1)/2;
    
    long unsigned stepper = 0;
    int r;
    
    if(trace==1){
      Rprintf("Using %i threads\n", number_of_threads);
      
      Rprintf("Progress:\n");
      for(r=0;r<101;r = r+5){
	if((r%10) == 0){ Rprintf("%d%%",r); }
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
    int t=0;
    while(t<nProfiles){
      strLine = DB(t);
      pProfile = new Profile(nLoci);
      for(nLoc=0;nLoc<nLoci;nLoc++){
	nPos = strLine.find(cDelim);
	strA1 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
	nPos = strLine.find(cDelim);
	strA2 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
	pProfile->m_ppdProfile[nLoc][0] = atoi(strA1.c_str());
	pProfile->m_ppdProfile[nLoc][1] = atoi(strA2.c_str());
      }
      vpProfiles.push_back(pProfile);
      t++;
    }
    
    int i, j;
    
    SEXP m;
    int m_size = ((1+wildcardEffect)*nLoci+1)*((1+wildcardEffect)*nLoci+1);
    PROTECT(m = allocVector(INTSXP, m_size));
    for(i=0;i<m_size;i++){ INTEGER(m)[i] = 0; }
    
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
    
    for(j=0; j<number_of_threads; j++){
      //Rprintf("Thread %i initiating... ", j);
      mthread_data = new thread_data();
      mthread_data->change_lock = &change_lock;
      mthread_data->i_queue_lock = &i_queue_lock;
      mthread_data->i = &i_row;
      mthread_data->id = j;
      mthread_data->nProfiles = nProfiles;
      mthread_data->nLoci = nLoci;
      mthread_data->hit = hit;
      mthread_data->trace = trace;
      mthread_data->single = single;
      mthread_data->wildcard = wildcard;
      mthread_data->wildcardEffect = wildcardEffect;
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
    
    if(nProfiles>=15 && trace==1){
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
    
    if(trace==1){
      for(j=progress_set; j<99; j++){
	Rprintf("=");
      }
    }
    
    for(j=0; j<number_of_threads; j++){
      rc = pthread_join(threads_container[j], (void **)&status);
      
      if (rc) {
	// REprintf("Error in thread joining: %i.\n", rc); 
	// fprintf(stderr, "Error in thread joining: %i.\n", rc);
	error("Error in thread joining"); 
	// exit(-1);
      }
    }
    
    if(trace==1){
      if(nProfiles<15){ 
	for(r=0;r<99;r++) Rprintf("=");
      }
      Rprintf("|\n");
    }
    
    //	for(i=0;i<(nLoci+1)*(nLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);
    
    SEXP r1, r2, mm, pp, fmm, fpp, list_names, list;
    PROTECT(r1 = allocVector(INTSXP, row1.size()));
    PROTECT(r2 = allocVector(INTSXP, row2.size()));
    PROTECT(mm = allocVector(INTSXP, match.size()));
    PROTECT(pp = allocVector(INTSXP, partial.size()));
    PROTECT(fmm = allocVector(INTSXP, fmatch.size()));
    PROTECT(fpp = allocVector(INTSXP, fpartial.size()));
    int matchlength = row1.size();
    for(i=0;i<matchlength;i++){
      INTEGER(r1)[i] = row1[i]; 
      INTEGER(r2)[i] = row2[i]; 
      INTEGER(mm)[i] = match[i]; 
      INTEGER(pp)[i] = partial[i];
      INTEGER(fmm)[i] = fmatch[i]; 
      INTEGER(fpp)[i] = fpartial[i];
    }
    
    PROTECT(list_names = allocVector(STRSXP,7));
    const char *names[7] = {"M","row1","row2","matches","partial","fmatches","fpartial"};
    for(i=0;i<7;i++) SET_STRING_ELT(list_names, i ,mkChar(names[i]));
    
    PROTECT(list = allocVector(VECSXP,7));
    SET_VECTOR_ELT(list, 0, m);
    SET_VECTOR_ELT(list, 1, r1);
    SET_VECTOR_ELT(list, 2, r2);
    SET_VECTOR_ELT(list, 3, mm);
    SET_VECTOR_ELT(list, 4, pp);
    SET_VECTOR_ELT(list, 5, fmm);
    SET_VECTOR_ELT(list, 6, fpp);
    setAttrib(list,R_NamesSymbol, list_names);
    
    UNPROTECT(9);
  
    return list;
#endif
  }
  
} // end extern
