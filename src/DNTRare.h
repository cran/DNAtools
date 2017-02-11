#include "Rcpp.h"
#include <map>
#include <vector>
#include <string>

using namespace Rcpp;
using namespace std;

class DNTRare {
  typedef double(DNTRare::*probFunc)(void);

  NumericVector m_vProbs; // original allele frequencies
  NumericVector m_vIsRare;
  int m_nAlleles;
  double m_dThreshold;
  double m_dTheta;
  
  // map that holds the probability function pointers for the combinations
  map<string, probFunc> m_mapProbFunc; 
  
  DNTRare(); // default constructor
  public:
  DNTRare(NumericVector q, NumericVector R, double r, double t); // explicit constructor
  
  private:  
    double pijklT(IntegerVector i);
    double pijkl(int *pnCounts, int *nCurr);
    double Pijkl(int i, int j, int k, int l);
    
    double pAAAR(void);
    double pAAAR_(void);
    double pAARA_(void);
    double pAARB_AB(void);
    double pAARB_BA(void);
    double pAABR_AB(void);
    double pAABR_BA(void);
    double pAARR(void);
    double pBARA(void);
    double pABRA(void);
    double pBAAR(void);
    double pABAR(void);
    double pABBR(void);
    double pBABR(void);
    double pABRB(void);
    double pBARB(void);
    double pABRC_ABC(void);
    double pABRC_ACB(void);
    double pABRC_CAB(void);
    double pABCR_ABC(void);
    double pABCR_ACB(void);
    double pABCR_CAB(void);
    double pABRR(void);
    double pARAR(void);
    double pARRA(void);
    double pRARA(void);
    double pARBR_AB(void);
    double pARBR_BA(void);
    double pARRB_AB(void);
    double pARRB_BA(void);
    double pRARB(void);
    double pARRR(void);
    double pRARR(void);
    double pRRRR(void);
    double pAAAA(void);
    double pAAAB(void);
    double pAABC(void);
    double pAABB(void);
    double pABAB(void);
    double pABAC(void);
    double pABCD(void);
    
    void setFunctionPointers(void);
    
  public:
    NumericVector prob(vector<string> vstrCombs);
};

