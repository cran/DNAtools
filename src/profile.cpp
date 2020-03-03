/*########################################################################*/
/* File: Profile.cpp                                                        */
/* Authors: James M. Curran, Torben Tvedebrink, and                       */
/*          Mikkel M. Andersen                                            */
/*                                                                        */
/* Version history                                                        */
/* ---------------------------------------------------------------------  */
/* Version  Date        Changes                          Author           */
/* -------  ----------  --------------------------       ---------------- */
/* 1.0      2013-10-23  Added first version number       JMC              */
/* 1.0      2013-10-23  Added new Profile constructor    JMC              */
/* 1.0      2013-10-23  Added Profile.toString method    JMC              */
/* 1.0      2013-10-23  Added extra arg and code         JMC              */  
/* 1.0                  to compare method to locus by                     */
/*                      scoring                                           */
/* 1.1      2018-03-22  Added vector constructor         MMA              */
/* 1.0-1    2013-10-24  Added wildcard scoring to        JMC              */
/*                      locus by locus scoring                            */
/* 1.0-1    2013-10-24  Fixed missing hom in Rare        JMC              */
/*                      swap (when B has R)                               */

#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>

// FIXME mikl
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#undef length

const int nRareCode = 990;
const int nWildcardCode = 0;

using namespace std;

class Profile{
public:
  int **m_ppnProfile;
  bool *m_pbRare;
  bool *m_pbWildcard;
  bool *m_pbHom;
  int m_nLoci;
  
public:
  // default constructor
  Profile(void){
  	m_ppnProfile = NULL;
  	m_pbRare = NULL;
  	m_pbWildcard = NULL;
  	m_pbHom = NULL;
  }

  // constructors 
  
  Profile(int nLoci){
  	m_nLoci = nLoci;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++)
  		m_ppnProfile[nLoc] = new int[2];
  }
  
  Profile(int *pnProfile, int nLoci){
  	m_nLoci = nLoci;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++){
  		m_ppnProfile[nLoc] = new int[2];
      m_ppnProfile[nLoc][0] = pnProfile[2 * nLoc];
      m_ppnProfile[nLoc][1] = pnProfile[2 * nLoc + 1];
      
      m_pbRare[nLoc] = (m_ppnProfile[nLoc][0] == nRareCode || m_ppnProfile[nLoc][1] == nRareCode) ? true : false;
  		m_pbWildcard[nLoc] = (m_ppnProfile[nLoc][0] == nWildcardCode || m_ppnProfile[nLoc][1] == nWildcardCode) ? true : false;
  		m_pbHom[nLoc] = (m_ppnProfile[nLoc][0] == m_ppnProfile[nLoc][1]) ? true : false;
  	}
  }

  // Does not check that profile.size() % 2 == 0
  Profile(const std::vector<int>& profile){
  	m_nLoci = profile.size() / 2;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++){
  		m_ppnProfile[nLoc] = new int[2];
      m_ppnProfile[nLoc][0] = profile[2 * nLoc];
      m_ppnProfile[nLoc][1] = profile[2 * nLoc + 1];
      
      m_pbRare[nLoc] = (m_ppnProfile[nLoc][0] == nRareCode || m_ppnProfile[nLoc][1] == nRareCode) ? true : false;
  		m_pbWildcard[nLoc] = (m_ppnProfile[nLoc][0] == nWildcardCode || m_ppnProfile[nLoc][1] == nWildcardCode) ? true : false;
  		m_pbHom[nLoc] = (m_ppnProfile[nLoc][0] == m_ppnProfile[nLoc][1]) ? true : false;
  	}
  }

  Profile(string &strLine, int nLoci, char cDelim = '\t'){
  	m_nLoci = nLoci;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	int nPos;
  	string strA1, strA2;
  
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++){
  		m_ppnProfile[nLoc] = new int[2];
  
  		nPos = strLine.find(cDelim);
  		strA1 = strLine.substr(0,nPos);
  		strLine = strLine.substr(nPos+1);
  		nPos = strLine.find(cDelim);
  		strA2 = strLine.substr(0,nPos);
  		strLine = strLine.substr(nPos+1);
  
  		m_ppnProfile[nLoc][0] = atoi(strA1.c_str());
  		m_ppnProfile[nLoc][1] = atoi(strA2.c_str());
  		m_pbRare[nLoc] = (m_ppnProfile[nLoc][0] == nRareCode || m_ppnProfile[nLoc][1] == nRareCode) ? true : false;
  		m_pbWildcard[nLoc] = (m_ppnProfile[nLoc][0] == nWildcardCode || m_ppnProfile[nLoc][1] == nWildcardCode) ? true : false;
  		m_pbHom[nLoc] = (m_ppnProfile[nLoc][0] == m_ppnProfile[nLoc][1]) ? true : false;
  	}
  }

  // copy constructor
  Profile(const Profile &p){
  	m_nLoci = p.m_nLoci;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++){
  		m_ppnProfile[nLoc] = new int[2];
  		m_ppnProfile[nLoc][0] = p.m_ppnProfile[nLoc][0];
  		m_ppnProfile[nLoc][1] = p.m_ppnProfile[nLoc][1];
  		m_pbRare[nLoc] = p.m_pbRare[nLoc];
  		m_pbWildcard[nLoc] = p.m_pbWildcard[nLoc];
  		m_pbHom[nLoc] = p.m_pbHom[nLoc];
  	}
  }

  // assignment operator overload
  const Profile& operator=(const Profile &p){
  	m_nLoci = p.m_nLoci;
  	m_ppnProfile = new int*[m_nLoci];
  	m_pbRare = new bool[m_nLoci];
  	m_pbWildcard = new bool[m_nLoci];
  	m_pbHom = new bool[m_nLoci];
  
  	int nLoc;
  	for(nLoc = 0; nLoc < m_nLoci; nLoc++){
  		m_ppnProfile[nLoc] = new int[2];
  		m_ppnProfile[nLoc][0] = p.m_ppnProfile[nLoc][0];
  		m_ppnProfile[nLoc][1] = p.m_ppnProfile[nLoc][1];
  		m_pbRare[nLoc] = p.m_pbRare[nLoc];
  		m_pbWildcard[nLoc] = p.m_pbWildcard[nLoc];
  		m_pbHom[nLoc] = p.m_pbHom[nLoc];
  	}
  
  	return *this;
  }

  void compare(Profile *pProf2, unsigned long &m2, unsigned long &m1, unsigned long &m0, 
                        unsigned long &fm2, unsigned long &fm1, 
    	                  bool bWildcardEffect = false, bool bRallele = false,
                        vector<int>* pvnScore = NULL){
  	int k;
  	int nA1, nA2, nB1, nB2, nR1, nR2;
  	bool bAHom, bBHom;
      
    for(k = 0; k < m_nLoci; k++){
      int nScore = 0;
      int nFScore = 0;
  
      nA1 = this->m_ppnProfile[k][0];
  		nA2 = this->m_ppnProfile[k][1];
  		bAHom = this->m_pbHom[k];
  
  		nB1 = pProf2->m_ppnProfile[k][0];
  		nB2 = pProf2->m_ppnProfile[k][1];
  		bBHom = pProf2->m_pbHom[k];
  
  		if(bRallele){
  			bool bAR = this->m_pbRare[k];
  			bool bBR = pProf2->m_pbRare[k];   
        
        //Rprintf("%d %d/%d %d/%d %c %c\n", k, nA1, nA2, nB1, nB2, (bAR ? 'T' : 'F'), (bBR ? 'T' : 'F'));
  
  			// Make comparions
  			if(!bAR && !bBR){
  				if(bAHom){ // Prof A hom
  					if(bBHom){ // Prof B hom
  						if(nA1 == nB1){ 
  							m2++; // AA,AA
                nScore = 2; 
  						}else{
  							m0++; // AA,BB
                nScore = 0;
  						}
  					}else{ // Prof B het
  						if(nA1 == nB1 || nA1 == nB2){ 
  							m1++; // AA,AB or BB,AB
                nScore = 1;
  						}else{ 
  							m0++;
                nScore = 0;
  						}
  					}
  				}else{ // Prof A het
  					if(bBHom){ // Prof B hom
  						if( nA1 == nB1 || nA2 == nB1){ //AB BB or BA AA
  							m1++;
                nScore = 1;
  						}else{
  							m0++;
                nScore = 0;
  						}
  					}else{ // 
  						if((nA1 == nB1 || nA1 == nB2) && (nA2 == nB1 || nA2 == nB2)){ 
  							m2++; // AB,AB or AB,BA
                nScore = 2;
  						}else if(((nA1 == nB1 || nA1 == nB2) && (nA2 != nB1 && nA2 != nB2)) || ((nA1 != nB1 && nA1 != nB2) && (nA2 == nB1 || nA2 == nB2))){
  							m1++; // AB,AC, AB,CA, AB,CB, ....
                nScore = 1;
  						}else{
  							m0++;
                nScore = 0;
  						}
  					}
  				}
  			}else if((bAR && !bBR) || (!bAR && bBR)){ /// RARE ALLELES IN ONLY ONE OF THE PROFILES
  				
  				bool bRHom = false;
  
  				if(bAR){ // Profile with rare is nR*
  					nR1 = nA1;
  					nR2 = nA2;
  
  					bRHom = bAHom;
  				}else{ // and other profile is b
  					nR1 = nB1;
  					nR2 = nB2;
  					nB1 = nA1;
  					nB2 = nA2;
  
  					bRHom = bBHom;
            bBHom = bAHom; // added v1.0-1 JMC
            
  //          if(bTrace)
  //            Rprintf("%d/%d %d/%d %c %c\n", nR1, nR2, nB1, nB2, (bRHom ? 'T' : 'F'), (bBHom ? 'T' : 'F'));
  				}
          
  
  				if(bRHom){
  					m2++; // Both rare variants - matches everything
            nScore = 2;
  				}else if(nR1 == nRareCode){  // The 'left' allele is rare
  					if(bBHom){ // Prof B hom
  						if(nR2 == nB1){ 
  							m1++; // RA,AA
                nScore = 1;
  						}else if(nR2 < nB1){ //****NOTE**** James: changed this to else if else. I don't think it will hurt but it may
  							m0++; // RA,BB
                nScore = 0;
  						}else{ // nR2 > nB1
  							m1++; // RB,AA
                nScore = 1;
  						}
  					}else{ // Prof B het
  						if(nR2 == nB1){ 
  							m1++; // RA,AB
                nScore = 1;
  						}else if(nR2 == nB2){
  							m2++; // RB,AB
                nScore = 2;
  						}else if(nR2 < nB1 && nR2 < nB2){
  							m0++; // RA,BC
                nScore = 0;
  						}else if(nR2 > nB1 || nR2 > nB2){
  							m1++; // RB,AC or RC,AB   
                nScore = 1;
  						}
  					}
  				}else if(nR2 == nRareCode){ // The 'right' allele is rare
  					if(bBHom){ // Prof B hom
  						if(nR1 == nB1){ // Changed this score to 2 v1.0-1 JMC and back to 1 0.1.18
  							m1++; // AR,AA
                nScore = 1;
  						}else if(nR1 < nB1){ //****NOTE**** James: changed this to else if else. I don't think it will hurt but it may
  							m1++; // AR,BB
                nScore = 1;
  						}else{ // nR1 > nB1
  							m0++; // BR,AA
                nScore = 0;
  						}
  					}else{ // Prof b het
  						if(nR1 == nB1){
  							m2++; // AR,AB
                nScore = 2;
  						}else if(nR1 == nB2){
  							m1++; // BR,AB
                nScore = 1;
  						}else if(nR1 < nB1 || nR1 < nB2){ // changed && to || 0.1-18
  							m1++; // AR,BC
                nScore = 1;
  						}else if(nR1 > nB1 || nR1 > nB2){
  							m0++; // BR,AC or CR,AB   
                nScore = 0;
  						}
  					}
  				}
  			}else if(bAR && bBR){ /// RARE ALLELES IN BOTH PROFILES
  				if(bAHom || bBHom){ // ((nA1 == nRareCode && nA2 == nRareCode) || (nB1 == nRareCode && nB2 == nRareCode)){ 
  					m2++; // a and/or b hom RR - matches everything
            nScore = 2;
  				}else if(nA2 == nRareCode){ // 'right' a allele rare AR:
  					if(nB2 == nRareCode){ // 'right' b allele rare: BR
  						if(nA1 == nB1){ 
  							m2++; // AR,AR
                nScore = 2;
  						}else{
  							m1++; // AR,BR
                nScore = 1;
  						}
  					}else{ // if(nB1 == nRareCode){ // 'left' b allele rare: RB //****NOTE**** James: changed this to else. I don't think it will hurt but it may
  						if(nA1 == nB2){
  							m1++; // AR,RA: Only partial matches because 7R' and R''7 would imply R'>7 and R''<7
                nScore = 1;
  						}else if(nA1 < nB2){ //****NOTE**** James: changed this to else if else. I don't think it will hurt but it may
  							m1++; // AR,BR: 7R' and 8R'' 
                nScore = 1;
  						}else{ // if(nA1 > nB2)
  							m0++; // BR,RA e.g. 9R' and R''7  means that R'>9 but R''<7
                nScore = 0;
  						}
  					}
  				}else if(nA1 == nRareCode){ // 'right' a allele rare RA:
  					if(nB2 == nRareCode){ // 'right' b allele rare: BR
  						if(nA2 == nB1){
  							m2++; // RA,AR
                nScore = 2;
  						}else if(nA2 < nB1){ //****NOTE**** James: changed this to else if else. I don't think it will hurt but it may
  							m0++; // RA,BR: 
                nScore = 0;
  						}else{ //if(nA2 > nB1)
  							m1++; // RB,AR
                nScore = 1;
  						}
  					}else{ // if(nB1 == nRareCode){ // 'left' b allele rare: RB
  						if(nA2 == nB2){
  							m2++; // AR,RA:
                nScore = 2;
  						}else{ 
  							m1++; // AR,RB
                nScore = 1;
  						}
  					}
  				}
  			}else{ 
  				//Rprintf("WHAT?\n"); 
  			} // Don't match cases: No R-alleles, R in one profiles, R in both profiles
  		}else{ // END R-allele // Otherwise
  			if(bAHom){ // Profile 1 is hom
  				if(nA1 == nWildcardCode){ // Both A alleles are wildcards "F"
  					fm2++;
            nFScore = 2;
  				}else{ // No wildcards in profile 1
  					if(bBHom){ // Profile 2 is hom
  						if(nB1 == nWildcardCode){
  							fm2++; // Both B alleles are wildcards "F"
  						}else if(nA1 == nB1){
  							m2++; // A genuine match
  							nScore = 2;
  							
                if(bWildcardEffect){
  								fm2++; // would be aF,aF
  								nScore = 2;
  							}
  						}else{
  							if(bWildcardEffect){
  								fm1++; // if analysing effect of F 
  								nScore = 1;
  							}
  						}
  					}else{ // Profile 2 is het
  						if(nB1 == nWildcardCode){ // Profile 2 is bF
  							if(nA1 == nB2){
  								fm2++;
  								nScore = 2;
  							}else{ 
  								fm1++;
  								nScore = 1;
  							}
  						}else{ // Profile 2 has no wildcards
  							if(nA1 == nB1 || nA1 == nB2){
  								m1++;
  								nScore = 1;
  								if(bWildcardEffect){
  									fm2++; // would be aF,ab
  									nFScore = 2;
  								}
  							}else{
  								if(bWildcardEffect)
  									fm1++;
  									nFScore = 1;
  							}
  						}
  					}
  				}
  			}else{ // Profile 1 is het
  				if(nA1 == nWildcardCode){ // Profile 1 is aF
  					if(bBHom){ // Profile 2 is hom
  						if(nB1 == nWildcardCode){
  							fm2++;
  							nFScore = 2;
  						}else if(nA2 == nB1){
  							fm2++;
  							nFScore = 2;
  						}else{
  							fm1++;
  							nFScore = 1;
  						}
  					}else{ // Profile 2 is het
  						if(nB1 == nWildcardCode){
  							if(nA2 == nB2){
  								fm2++;
  								nFScore = 2;
  							}else{
  								fm1++;
  								nFScore = 1;
  							}
  						}else{ // Profile 2 has no wildcards
  							if(nA2 == nB1 || nA2 == nB2){ 
  								fm2++;
  								nFScore = 2;
  							}else{
  								fm1++;
  								nFScore = 1;
  							}
  						}
  					}
  				}else{ // profile 1 has no wildcards
  					if(bBHom){ // Profile 2 is hom
  						if(nB1 == nWildcardCode){
  							fm2++;
  							nFScore = 2;
  						}else if(nA1 == nB1 || nA2 == nB1){ 
  							m1++;
  							nScore = 1;
  							
                if(bWildcardEffect){
  								fm2++;
  								nFScore = 2;
  							}
  						}else{
  							if(bWildcardEffect){
  								fm1++;
  								nFScore = 1;
  							}
  						}
  					}else{ // Profile 2 is het
  						if(nB1 == nWildcardCode){
  							if(nA1 == nB2 || nA2 == nB2){
  								fm2++;
  								nFScore = 2;
  							}else{
  								fm1++;
  								nFScore = 1;
  							}
  						}else{ // None of the profiles has wildcards
  							if(nA1 == nB1 && nA2 == nB2){
  								m2++;
  								nScore = 2;
  								if(bWildcardEffect){
  									fm2++;
  									nFScore = 2;
  								}
  							}else if((nA1 == nB1 && nA2 != nB2) || (nA1 == nB2 && nA2 != nB1) || (nA1 != nB1 && nA2 == nB2) || (nA1 != nB2 && nA2 == nB1)){
  								m1++;
  								nScore = 1;
  								
                  if(bWildcardEffect){
  									fm1++;
  									nFScore = 1;
  								}
  							}
  						}
  					}
  				}
  			}
  		} // End else (otherwise)
      if(pvnScore != NULL){
        (*pvnScore)[k] = bWildcardEffect ? nFScore : nScore;  // changed to deal with F designations 2013-10-24 JMC
      }
  	} // End k
    
    //return result;
  }
    
  string toString(bool bShowHoms){
    ostringstream oss;
    
    for(int k = 0; k < m_nLoci; k++){
      oss << this->m_ppnProfile[k][0] << '/' << this->m_ppnProfile[k][1];
      
      if(bShowHoms)
        oss << ' ' << (m_pbHom[k] ? 'T' : 'F');
        
      if(k < (m_nLoci - 1))
        oss << ' ';
    }
    
    return oss.str();
  }
};
