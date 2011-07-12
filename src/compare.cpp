#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

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

  int unsigned comps = nProfiles*(nProfiles-1)/2;
  if(single==1) comps = (nProfiles-1);

  int unsigned stepper = 0;
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

} // end extern
