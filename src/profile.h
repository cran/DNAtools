/*########################################################################*/
/* File: Profile.h                                                        */
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
/* 1.0      2013-10-23  Added extra arg to compare       JMC              */  
/* 1.0                  method for locus by locus                         */
/*                      scoring                                           */
/* 1.1      2018-03-22  Added vector constructor         MMA              */

#include <string>
#include <sstream>
#include <vector>

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
	Profile(void);
  
  // constructors
	Profile(int nLoci);
  Profile(int *pnProfile, int nLoci);
  Profile(const std::vector<int>& profile);
	Profile(string &strLine, int nLoci, char cDelim = '\t');

  // copy constructor
  Profile(const Profile &p);
  
  // assignment operator overload
	const Profile& operator=(const Profile &p);

  // compare two profiles
	void compare(Profile *pProf2, unsigned long &m2, unsigned long &m1, unsigned long &m0, unsigned long &fm2, unsigned long &fm1, 
		         bool bWildcardEffect = false, bool bRallele = false, vector<int>* pvnScore = NULL);
  
  // printing code
  string toString(bool bShowHoms = false);
};

