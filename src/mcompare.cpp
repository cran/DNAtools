
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <queue>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>


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
#ifndef WIN32
	usleep(500000);
#endif
}

extern "C" { // begin extern

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
	string m_strName;
	double **m_ppdProfile;
	int m_nLoci;

public:
	Profile(void){
		m_strName.clear();
		m_ppdProfile = NULL;
	};

	Profile(string &strName, int nLoci){
		m_strName = strName;
		m_nLoci = nLoci;

		m_ppdProfile = new double*[nLoci];

		int nLoc;
		for(nLoc=0;nLoc<nLoci;nLoc++)
			m_ppdProfile[nLoc] = new double[2];
	};

	Profile(const Profile &p){
		m_strName = p.m_strName;
		m_nLoci = p.m_nLoci;

		m_ppdProfile = new double*[m_nLoci];

		int nLoc;
		for(nLoc=0;nLoc<m_nLoci;nLoc++){
			m_ppdProfile[nLoc] = new double[2];
			m_ppdProfile[nLoc][0] = p.m_ppdProfile[nLoc][0];
			m_ppdProfile[nLoc][1] = p.m_ppdProfile[nLoc][1];
		}
	};

	const Profile& operator=(const Profile &p){
		m_strName = p.m_strName;
		m_nLoci = p.m_nLoci;

		m_ppdProfile = new double*[m_nLoci];

		int nLoc;
		for(nLoc=0;nLoc<m_nLoci;nLoc++){
			m_ppdProfile[nLoc] = new double[2];
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
	int nLoci;
	int hit;
	int trace;
	int single;
	long unsigned* stepper;
	vector<Profile*> *vpProfiles;
	SEXP m;
	vector<int> *row1;
	vector<int> *row2;
	vector<int> *match;
	vector<int> *partial;
};
#endif

#ifdef DBCOMPARE_MULTICORE
void *compare_profiles(void *threadarg) {
	thread_data *my_data;
	my_data = (thread_data*) threadarg;

	pthread_mutex_t* change_lock = my_data->change_lock;
	pthread_mutex_t* i_queue_lock = my_data->i_queue_lock;

	int id = my_data->id;
	int nProfiles = my_data->nProfiles;
	int nLoci = my_data->nLoci;
	int hit = my_data->hit;
	int trace = my_data->trace;
	int single = my_data->single;
	long unsigned *stepper = my_data->stepper;
	vector<Profile*> *vpProfiles = my_data->vpProfiles;
	SEXP m = my_data->m;
	vector<int> *row1 = my_data->row1;
	vector<int> *row2 = my_data->row2;
	vector<int> *match = my_data->match;
	vector<int> *partial = my_data->partial;

	int i;

	vector<int> local_row1;
	vector<int> local_row2;
	vector<int> local_match;
	vector<int> local_partial;

	int j, k, m1, m2;
	double dA1, dA2, dB1, dB2;
	
	int m_size = (nLoci+1)*(nLoci+1);
	
	SEXP local_m;
  PROTECT(local_m = allocVector(INTSXP, m_size));
	
	for(j=0;j<m_size;j++){ INTEGER(local_m)[j] = 0; }

	Profile *pProf1, *pProf2;

	for (;;) { // we break when there's no more i's left
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
			// if a single comparison is needed:
			if(single==1) i = nProfiles;

			pProf2 = (*vpProfiles)[j];

			m2 = 0;
			m1 = 0;

			for(k=0;k<nLoci;k++){
				dA1 = pProf1->m_ppdProfile[k][0];
				dA2 = pProf1->m_ppdProfile[k][1];
				dB1 = pProf2->m_ppdProfile[k][0];
				dB2 = pProf2->m_ppdProfile[k][1];

				if(dA1!=0 && dA2!=0 && dB1!=0 && dB2!=0){
					if(dA1==dB1 && dA2 == dB2){
						m2++;
					}else if((dA1==dB1 && dA2!=dB2) || (dA1==dB2 && dA2!=dB1) || (dA1!=dB1 && dA2==dB2) || (dA1!=dB2 && dA2==dB1)){
						m1++;
					}
				}
			}

			// thread_m(m2,m1)++;
			INTEGER(local_m)[m2*(nLoci+1) + m1] ++;			

			if(m2>=hit){
				//Rprintf("hit = %i, m2 = %i\n", hit, m2);
				//	 prof1.push_back(pProf1->m_strName);
				//	 prof2.push_back(pProf2->m_strName);
				if(single==1) local_row1.push_back(1); // Always first row when making comparisons to a given profile
				else local_row1.push_back(i+1);
				local_row2.push_back(j+1);
				local_match.push_back(m2);
				local_partial.push_back(m1);
			}

		} // end for(j)

		pthread_mutex_lock(change_lock); // lock because we need to write to the shared data
		*stepper += nProfiles - (i+1);
		pthread_mutex_unlock(change_lock); // unlock again
	}

	int n1 = local_row1.size();

	pthread_mutex_lock(change_lock); // lock because we need to write to the shared data
	for(j=0;j<n1;j++) {
		row1->push_back(local_row1[j]);
		row2->push_back(local_row2[j]);
		match->push_back(local_match[j]);
		partial->push_back(local_partial[j]);
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

  vector<Profile*> vpProfiles;
  int nProfiles = DB.size();

  string strLine;

  int nLoc;
  string strID, strA1, strA2;
  size_t nPos;
  Profile *pProfile;
  char cDelim = '\t';

  long unsigned comps = nProfiles*(nProfiles-1)/2;
  if(single==1) comps = (nProfiles-1);

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

  int t=0;
  while(t<nProfiles){
    strLine = DB(t);
    nPos = strLine.find(cDelim);
    strID = strLine.substr(0,nPos);
    strLine = strLine.substr(nPos+1);
    
    pProfile = new Profile(strID, nLoci);

    for(nLoc=0;nLoc<nLoci;nLoc++){
      nPos = strLine.find(cDelim);
      strA1 = strLine.substr(0,nPos);
      strLine = strLine.substr(nPos+1);
      
//       if(nLoc<(nLoc-1)){				
// 	strA2 = strLine;
//       }else{
	nPos = strLine.find(cDelim);
	strA2 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
//       }
      
      if(atof(strA1.c_str())>atof(strA2.c_str())){
	string strTmp = strA1;
	strA1 = strA2;
	strA2 = strTmp;
      }
      
      pProfile->m_ppdProfile[nLoc][0] = atof(strA1.c_str());
      pProfile->m_ppdProfile[nLoc][1] = atof(strA2.c_str());
    }
    
    vpProfiles.push_back(pProfile);
    t++;
  }

  int i, j, k, m1, m2;
  
  SEXP m;
  PROTECT(m = allocVector(INTSXP, (nLoci+1)*(nLoci+1)));
  for(i=0;i<(nLoci+1)*(nLoci+1);i++){ INTEGER(m)[i] = 0; }

//   vector<string> prof1;
//   vector<string> prof2;
  vector<int> row1;
  vector<int> row2;
  vector<int> match;
  vector<int> partial;

  double dA1, dA2, dB1, dB2;
	
  Profile *pProf1, *pProf2;

  for(i=0;i<nProfiles;i++){
    pProf1 = vpProfiles[i];
    for(j=i+1;j<nProfiles;j++){
      // if a single comparison is needed:
      if(single==1) i = nProfiles;

      pProf2 = vpProfiles[j];
	    
      m2 = 0;
      m1 = 0;
	    
      for(k=0;k<nLoci;k++){
	dA1 = pProf1->m_ppdProfile[k][0];
	dA2 = pProf1->m_ppdProfile[k][1];
	dB1 = pProf2->m_ppdProfile[k][0];
	dB2 = pProf2->m_ppdProfile[k][1];
	
	if(dA1!=0 && dA2!=0 && dB1!=0 && dB2!=0){
	  if(dA1==dB1 && dA2 == dB2){
	    m2++;
	  }else if((dA1==dB1 && dA2!=dB2) || (dA1==dB2 && dA2!=dB1) || (dA1!=dB1 && dA2==dB2) || (dA1!=dB2 && dA2==dB1)){
	    m1++;
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
      INTEGER(m)[m2*(nLoci+1)+m1] ++;
      
      if(m2>=hit){
// 	prof1.push_back(pProf1->m_strName);
// 	prof2.push_back(pProf2->m_strName);
	if(single==1) row1.push_back(1); // Always first row when making comparisons to a given profile
	else row1.push_back(i+1);
	row2.push_back(j+1);
	match.push_back(m2);
	partial.push_back(m1);
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
  
  SEXP r1, r2, mm, pp, list_names, list;
  PROTECT(r1 = allocVector(INTSXP, row1.size()));
  PROTECT(r2 = allocVector(INTSXP, row2.size()));
  PROTECT(mm = allocVector(INTSXP, match.size()));
  PROTECT(pp = allocVector(INTSXP, partial.size()));
  int matchlength = row1.size();
  for(i=0;i<matchlength;i++){
    INTEGER(r1)[i] = row1[i]; 
    INTEGER(r2)[i] = row2[i]; 
    INTEGER(mm)[i] = match[i]; 
    INTEGER(pp)[i] = partial[i]; 
  }

  PROTECT(list_names = allocVector(STRSXP,5));
  const char *names[5] = {"M","row1","row2","matches","partial"};
  for(i=0;i<5;i++) SET_STRING_ELT(list_names, i ,mkChar(names[i]));

  PROTECT(list = allocVector(VECSXP,5));
  SET_VECTOR_ELT(list, 0, m);
  SET_VECTOR_ELT(list, 1, r1);
  SET_VECTOR_ELT(list, 2, r2);
  SET_VECTOR_ELT(list, 3, mm);
  SET_VECTOR_ELT(list, 4, pp);
  setAttrib(list,R_NamesSymbol, list_names);
  
  UNPROTECT(7);
  
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
	if(single==1) comps = (nProfiles-1);

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

	int t=0;
	while(t<nProfiles){
		strLine = DB(t);
		nPos = strLine.find(cDelim);
		strID = strLine.substr(0,nPos);
		strLine = strLine.substr(nPos+1);

		pProfile = new Profile(strID, nLoci);

		for(nLoc=0;nLoc<nLoci;nLoc++){
			nPos = strLine.find(cDelim);
			strA1 = strLine.substr(0,nPos);
			strLine = strLine.substr(nPos+1);

//			 if(nLoc<(nLoc-1)){
//	 strA2 = strLine;
//			 }else{
	nPos = strLine.find(cDelim);
	strA2 = strLine.substr(0,nPos);
	strLine = strLine.substr(nPos+1);
//			 }

			if(atof(strA1.c_str())>atof(strA2.c_str())){
	string strTmp = strA1;
	strA1 = strA2;
	strA2 = strTmp;
			}

			pProfile->m_ppdProfile[nLoc][0] = atof(strA1.c_str());
			pProfile->m_ppdProfile[nLoc][1] = atof(strA2.c_str());
		}

		vpProfiles.push_back(pProfile);
		t++;
	}

	int i, j, m1, m2;

	SEXP m;
	PROTECT(m = allocVector(INTSXP, (nLoci+1)*(nLoci+1)));
	for(i=0;i<(nLoci+1)*(nLoci+1);i++){ INTEGER(m)[i] = 0; }

//	 vector<string> prof1;
//	 vector<string> prof2;
	vector<int> row1;
	vector<int> row2;
	vector<int> match;
	vector<int> partial;

	double dA1, dA2, dB1, dB2;

	Profile *pProf2;

	vector<thread_data> thread_data_array;

	int rc, status;
	pthread_mutex_t change_lock;
	pthread_mutex_t i_queue_lock;

	if (pthread_mutex_init(&change_lock, NULL)) {
		fprintf(stderr, "Could not initialize change_lock mutex, aborting.\n");
		exit(1);
	}

	if (pthread_mutex_init(&i_queue_lock, NULL)) {
		fprintf(stderr, "Could not initialize i_queue_lock mutex, aborting.\n");
		exit(1);
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
		mthread_data->stepper = &stepper;
		mthread_data->vpProfiles = &vpProfiles;
		mthread_data->m = m;
		mthread_data->row1 = &row1;
		mthread_data->row2 = &row2;
		mthread_data->match = &match;
		mthread_data->partial = &partial;
		thread_data_array.push_back(*mthread_data);
		
		pthread_t thr;
		rc = pthread_create(&thr, NULL, compare_profiles, (void *) mthread_data);
		threads_container.push_back(thr);

		if (rc) {
			fprintf(stderr, "Error in thread creation: %i.\n", rc);
			exit(-1);
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
	for(j=progress_set; j<99; j++){
		Rprintf("=");
	}

	for(j=0; j<number_of_threads; j++){
		rc = pthread_join(threads_container[j], (void **)&status);

		if (rc) {
			fprintf(stderr, "Error in thread joining: %i.\n", rc);
			exit(-1);
		}
	}

	if(trace==1){
		if(nProfiles<15){ 
			for(r=0;r<99;r++) Rprintf("=");
		}
		Rprintf("|\n");
	}

	//	for(i=0;i<(nLoci+1)*(nLoci+2)/2;i++) Rprintf("%d ",INTEGER(m)[i]);

	SEXP r1, r2, mm, pp, list_names, list;
	PROTECT(r1 = allocVector(INTSXP, row1.size()));
	PROTECT(r2 = allocVector(INTSXP, row2.size()));
	PROTECT(mm = allocVector(INTSXP, match.size()));
	PROTECT(pp = allocVector(INTSXP, partial.size()));
	int matchlength = row1.size();
	for(i=0;i<matchlength;i++){
		INTEGER(r1)[i] = row1[i]; 
		INTEGER(r2)[i] = row2[i]; 
		INTEGER(mm)[i] = match[i]; 
		INTEGER(pp)[i] = partial[i];
	}

	PROTECT(list_names = allocVector(STRSXP,5));
	const char *names[5] = {"M","row1","row2","matches","partial"};
	for(i=0;i<5;i++) SET_STRING_ELT(list_names, i ,mkChar(names[i]));

	PROTECT(list = allocVector(VECSXP,5));
	SET_VECTOR_ELT(list, 0, m);
	SET_VECTOR_ELT(list, 1, r1);
	SET_VECTOR_ELT(list, 2, r2);
	SET_VECTOR_ELT(list, 3, mm);
	SET_VECTOR_ELT(list, 4, pp);
	setAttrib(list,R_NamesSymbol, list_names);

	UNPROTECT(7);

	return list;
	#endif
}

} // end extern
