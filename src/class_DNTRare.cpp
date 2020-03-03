#ifndef _CLASS_DNTRARE
#define _CLASS_DNTRARE 1

#include <Rcpp.h>

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
  
  // default constructor
  DNTRare() {
    m_dThreshold = 0;
    m_dTheta = 0;
    m_nAlleles = 0;
  }
  
public:
  // explicit constructor
  DNTRare(NumericVector q, NumericVector R, double r, double t) {
    
    setFunctionPointers();
    
    m_dThreshold = r;
    m_dTheta = t;
    m_nAlleles = q.size();
    
    // re-scale q
    double dSum = sum(q);
    int i;
    
    for(i = 0; i < m_nAlleles; i++)
      q[i] /= dSum;
    
    m_vProbs = NumericVector(q);
    
    // Inserting rare allele probabilities (R- and R+)
    if(R.size() ==1) {
      m_vProbs.push_front(R[0] * 0.5);
      m_vProbs.push_back(R[0] * 0.5);
    } else {
      m_vProbs.push_front(R[0]);
      m_vProbs.push_back(R[1]);
    }
    
    double RR = sum(R);
    m_vIsRare = LogicalVector(m_nAlleles, false);
    m_vIsRare.push_front(true); // First and last are R- and R+
    m_vIsRare.push_back(true);
    
    // Scaling probabilitity vector and identifies rares
    for(i = 1; i <= m_nAlleles; i++) {
      m_vProbs[i] *= 1 - RR; //= q[i-1]*(1-RR);
      m_vIsRare[i] = (q[i - 1] < m_dThreshold) ? true : false;
    }
  }
  
private:
  double pijklT(IntegerVector i){
    int n = i.size();
    double res = 0.0;
    if(n>1){
      IntegerVector j(n-1);
      int x = 0;
      for(int k=1; k<n; k++){
        j[k-1] = i[k];
        if(i[0]==i[k]) x++;
      }
      
      Rprintf("j: ");
      for(int k = 0; k < n - 1; k++){
        Rprintf("%d ", j[k]);
      }
      Rprintf("\ni: ");
      for(int k = 0; k < n; k++){
        Rprintf("%d ", i[k]);
      }
      Rprintf("\nx %d\n", x);
      res = pijklT(j)*(x * m_dTheta + (1 - m_dTheta)* m_vProbs[i[0]]);
    }
    else res = m_vProbs[i[0]];
    return res;
  }
  
  
  double pijkl(int *pnCounts, int *nCurr) {
    int m = 4 - (*nCurr);
    double dResult = 0;
    
    if(m > 1) {
      int *pnFirst  = pnCounts + 1; //&(pnCounts[1]);
      *nCurr += 1;
      int x = 0;
      
      for(int k = 1; k < m; k++) {
        if(pnCounts[0] == pnCounts[k]) {
          x++;
        }
      }
      //         Rprintf("j: ");
      //    for(int k = 0; k < m - 1; k++){
      //      Rprintf("%d ", pnFirst[k]);
      //    }
      //    Rprintf("\ni: ");
      //    for(int k = 0; k < m; k++){
      //      Rprintf("%d ", pnCounts[k]);
      //    }
      //        Rprintf("x %d\n", x);
      dResult = pijkl(pnFirst, nCurr) * ( x * m_dTheta + (1 - m_dTheta) * m_vProbs[pnCounts[0]]);
    } else {
      dResult = m_vProbs[pnCounts[0]];
    }
    
    return dResult;
    /*int n = i.size();
    double res = 0.0;
    if(n>1){
    IntegerVector j(n-1);
    int x = 0;
    for(int k = 1; k < n; k++){
    j[k-1] = i[k];
    if(i[0]==i[k]) x++;
    }
    res = pijkl(j,p,t)*(x*t+(1-t)*p[i[0]]);
    }
    else res = p[i[0]];
    return res;*/
  }
  
  double Pijkl(int i, int j, int k, int l) {
    double D = (1 + 2 * m_dTheta) * (1 + m_dTheta);
    //IntegerVector y = IntegerVector::create(i, j, k, l);
    int nCounts[4] = {i, j, k, l};
    int nCurr = 0;
    
    //Rprintf("Combination: %d %d %d %d\n", i, j, k, l);
    
    //double dResult = pijklT(y);
    double dResult = pijkl(&nCounts[0], &nCurr);
    
    return dResult / D;
  }
  
  // probability function definitions
  double pAAAR(void) { // match 2
#ifdef _DEBUG
    Rprintf("pAAAR\n");
#endif
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 0; j <= m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            dResult += Pijkl(i,i,i,j) + Pijkl(i,i,j,i) + Pijkl(j,i,i,i) + Pijkl(i,j,i,i);
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAAAR_(void) { // match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            dResult += 2*Pijkl(i,i,i,j);
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAARA_(void) { // match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 0; j < i; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            dResult += 2*Pijkl(i,i,i,j);
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  
  double pAARB_AB(void) { // A<B & R<B match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = 0; k < j; k++) {
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,i,k,j) + Pijkl(i,i,j,k);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAARB_BA(void) { // B<A & R<B match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = 0; k < j; k++) {
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,i,k,j) + Pijkl(i,i,j,k);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAABR_AB(void) { // A<B & B<R match 0
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = j + 1; k <=  m_nAlleles + 1; k++) {
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,i,k,j) + Pijkl(k,j,i,i);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAABR_BA(void) { // B<A & B<R match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = j + 1; k <=  m_nAlleles + 1; k++) {
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,i,k,j) + Pijkl(k,j,i,i);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAARR(void) { // B<A & B<R match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 0; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) {
            for(int k = 0; k <=  m_nAlleles + 1; k++) {
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,i,k,j)  + Pijkl(k,j,i,i);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  double pBARA(void) { // B<A, R<A match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // B = A_j with B<A
          if(!m_vIsRare[j]) {
            for(int k = 0; k < i; k++) { // R = A_k with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += 3*(Pijkl(j,i,k,i) + Pijkl(k,i,j,i));
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABRA(void) { // A<B, R<A match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // B = A_j with B>A
          if(!m_vIsRare[j]) {
            for(int k = 0; k < i; k++) { // R = A_k with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += 3*(Pijkl(i,j,k,i) + Pijkl(k,i,i,j));
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pBAAR(void) { // B<A, A<R match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // B = A_j with B<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_k with R>A
              if(m_vIsRare[k]) { // A_k is a rare allele
                /* if(p[k]>0){
                std::cout<<i<<" "<<j<<" "<<i<<" "<<k<<"\n";
                std::cout<<i<<" "<<j<<" "<<k<<" "<<i<<"\n";
                std::cout<<i<<" "<<k<<" "<<j<<" "<<i<<"\n";
                std::cout<<i<<" "<<k<<" "<<i<<" "<<j<<"\n";
                std::cout<<k<<" "<<i<<" "<<i<<" "<<j<<"\n";
                std::cout<<k<<" "<<i<<" "<<j<<" "<<i<<"\n";
              } */
                dResult += Pijkl(j,i,i,k)*6;// + Pijkl(i,k,j,i);
            }
          }
        }
      }
    }
  }
    return dResult;
  }
  
  
  double pABAR(void) { // A<B, A<R match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // B = A_j with B<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_k with R>A
              if(m_vIsRare[k]) { // A_k is a rare allele
                /*	    if(p[k]>0){
                std::cout<<i<<" "<<j<<" "<<i<<" "<<k<<"\n";
                std::cout<<i<<" "<<j<<" "<<k<<" "<<i<<"\n";
                //	      std::cout<<j<<" "<<i<<" "<<i<<" "<<k<<"\n";
                //	      std::cout<<j<<" "<<i<<" "<<k<<" "<<i<<"\n";
                //
                std::cout<<i<<" "<<k<<" "<<i<<" "<<j<<"\n";
                std::cout<<k<<" "<<i<<" "<<i<<" "<<j<<"\n";
                std::cout<<i<<" "<<k<<" "<<j<<" "<<i<<"\n";
                std::cout<<k<<" "<<i<<" "<<j<<" "<<i<<"\n";
              }*/
                dResult += Pijkl(i,j,i,k)*6;
            }
          }
        }
      }
    }
  }
    return dResult;
}
  
  
  double pABBR(void) { // A<B, B<R match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(!m_vIsRare[j]) { // A_i!=A_j
            for(int k = j + 1; k <=  m_nAlleles + 1; k++) {
              if(m_vIsRare[k]) {
                dResult += Pijkl(i,j,j,k) + Pijkl(j,k,i,j);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pBABR(void) { // B<A, B<R match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // B = A_j with B!=A
          if(!m_vIsRare[j]) { // A_i!=A_j
            for(int k = j + 1; k <=  m_nAlleles + 1; k++) {
              if(m_vIsRare[k]) {
                dResult += Pijkl(i,j,j,k) + Pijkl(j,k,i,j);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABRB(void) { // A<B, B<R match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(!m_vIsRare[j]) { // A_i!=A_j
            for(int k = 0; k < j; k++) {
              if(m_vIsRare[k]) {
                dResult += Pijkl(i,j,j,k) + Pijkl(j,k,i,j);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pBARB(void) { // B<A, B<R match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // B = A_j with B!=A
          if(!m_vIsRare[j]) { // A_i!=A_j
            for(int k = 0; k < j; k++) {
              if(m_vIsRare[k]) {
                dResult += Pijkl(i,j,j,k) + Pijkl(j,k,i,j);
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  double pABRC_ABC(void) { // A!=B, R<C, {A,B}<C match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::max(i,j)+1; k <= m_nAlleles; k++) { // C = A_k with {A,B}<C
              if(!m_vIsRare[k]) {
                for(int l = 0; l < k; l++) { // R = A_l with R<C
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,l,k) + Pijkl(l,k,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABRC_ACB(void) { // A!=B, R<C, {A<C<B or B<C<A} match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::min(i,j)+1; k < std::max(i,j); k++) { // C = A_k with A < C < B or B < C < A
              if(!m_vIsRare[k]) {
                for(int l = 0; l < k; l++) { // R = A_l with R<C
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,l,k) + Pijkl(l,k,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABRC_CAB(void) { // A!=B, R<C, C<{A,B} match 0
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::min(std::min(1,j),std::min(1,i)); k < std::min(i,j); k++) { // C = A_k with C < A < B or C < B < A
              if(!m_vIsRare[k]) {
                for(int l = 0; l < k; l++) { // R = A_l with R<C
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,l,k) + Pijkl(l,k,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  double pABCR_ABC(void) { // A!=B, C<R, {A,B}<C match 0
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::max(i,j)+1; k <= m_nAlleles; k++) { // C = A_k with {A,B}<C
              if(!m_vIsRare[k]) {
                for(int l = k + 1; l <=  m_nAlleles + 1; l++) { // R = A_l with C<R
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,k,l) + Pijkl(k,l,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABCR_ACB(void) { // A!=B, C<R, {A<C<B or B<C<A} match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::min(i,j)+1; k < std::max(i,j); k++) { // C = A_k with A < C < B or B < C < A
              if(!m_vIsRare[k]) {
                for(int l = k + 1; l <=  m_nAlleles + 1; l++) { // R = A_l with C<R
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,k,l) + Pijkl(k,l,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABCR_CAB(void) { // A!=B, C<R, C<{A,B} match 1
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = std::min(std::min(1,j),std::min(1,i)); k < std::min(i,j); k++) { // C = A_k with C < A < B or C < B < A
              if(!m_vIsRare[k]) {
                for(int l = k + 1; l <=  m_nAlleles + 1; l++) { // R = A_l with R<C
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,k,l) + Pijkl(k,l,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABRR(void) { // A!=B, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // B = A_j with B!=A
          if(i!=j && !m_vIsRare[j]) { // A_i!=A_j
            for(int k = 0; k <=  m_nAlleles + 1; k++) { // R = A_k
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l <=  m_nAlleles + 1; l++) { // R = A_l
                  if(m_vIsRare[l]) { // A_l is a rare allele
                    dResult += Pijkl(i,j,k,l) + Pijkl(k,l,i,j);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARAR(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,j,i,k)*4;
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARRA(void) { // R<A, match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = 0; k < i; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += (Pijkl(i,j,i,k) + Pijkl(i,k,i,j))*4;
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pRARA(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 0; j < i; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = 0; k < i; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                dResult += Pijkl(i,j,i,k)*4;
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  double pARBR_AB(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = j + 1; l <=  m_nAlleles + 1; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,j,l) + Pijkl(j,l,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARBR_BA(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = j + 1; l <=  m_nAlleles + 1; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,j,l) + Pijkl(j,l,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARRB_AB(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l < j; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,l,j) + Pijkl(l,j,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARRB_BA(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j < i; j++) { // R = A_j with R<A
          if(!m_vIsRare[j]) {
            for(int k = i + 1; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l < j; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,j,l) + Pijkl(j,l,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pRARB(void) { // R<A, match 2
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 1; j <= m_nAlleles; j++) { // R = A_j with R<A
          if(i!=j && !m_vIsRare[j]) {
            for(int k = 0; k < i; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l < j; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += 2*(Pijkl(i,k,j,l) + Pijkl(j,l,i,k));
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pARRR(void) { // R<A, match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = i + 1; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = 0; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l <=  m_nAlleles + 1; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,j,l) + Pijkl(j,l,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pRARR(void) { // R<A, match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) {
        for(int j = 0; j <= i; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = 0; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l <=  m_nAlleles + 1; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += (Pijkl(i,k,j,l) + Pijkl(j,l,i,k))*2;
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  
  double pRRRR(void) { // R<A, match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 0; i <=  m_nAlleles + 1; i++) { // A = A_i
      if(m_vIsRare[i]) { // A_k is a rare allele
        for(int j = 0; j <=  m_nAlleles + 1; j++) { // R = A_j with R<A
          if(m_vIsRare[j]) { // A_k is a rare allele
            for(int k = 0; k <=  m_nAlleles + 1; k++) { // R = A_j with R<A
              if(m_vIsRare[k]) { // A_k is a rare allele
                for(int l = 0; l <=  m_nAlleles + 1; l++) {
                  if(m_vIsRare[l]) { // A_k is a rare allele
                    dResult += Pijkl(i,k,j,l);
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  //OK
  
  double pAAAA(void) { // R<A, match 2
    
    double dResult = 0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i]) dResult += Pijkl(i,i,i,i);
    }
    return dResult;
  }
  
  
  double pAAAB(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) dResult += Pijkl(i,i,i,j) + Pijkl(i,i,j,i) + Pijkl(i,j,i,i) + Pijkl(j,i,i,i);
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAABC(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) {
              for(int k = 1; k <= m_nAlleles; k++) {
                if(!m_vIsRare[k])
                {
                  if(i!=k && j!=k) dResult += Pijkl(i,i,j,k) + Pijkl(j,k,i,i);
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pAABB(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) dResult += Pijkl(i,i,j,j);
          }
        }
      }
    }
    return dResult;
  }
  
  
  // OK
  
  double pABAB(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) {
              dResult += Pijkl(i,j,i,j) + Pijkl(i,j,j,i);
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABAC(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) {
              for(int k = 1; k <= m_nAlleles; k++) {
                if(!m_vIsRare[k])
                {
                  if(i!=k && j!=k) dResult += Pijkl(i,j,i,k) + Pijkl(i,j,k,i) + Pijkl(j,i,k,i) + Pijkl(j,i,i,k);
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  double pABCD(void) { // R<A, match 2
    
    double dResult = 0.0;
    // Note: p = (R-,q1,q2,...,qn,R+)
    for(int i = 1; i <= m_nAlleles; i++) { // A = A_i
      if(!m_vIsRare[i])
      {
        for(int j = 1; j <= m_nAlleles; j++) {
          if(!m_vIsRare[j])
          {
            if(i!=j) {
              for(int k = 1; k <= m_nAlleles; k++) {
                if(!m_vIsRare[k])
                {
                  if(i!=k && j!=k) {
                    for(int l = 1; l <= m_nAlleles; l++) {
                      if(!m_vIsRare[l])
                      {
                        if(i!=l && j!=l && k!=l) dResult += Pijkl(i,j,k,l);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return dResult;
  }
  
  
  void setFunctionPointers(void){
    m_mapProbFunc["AAAR"] = &DNTRare::pAAAR;
    m_mapProbFunc["AAAR_"] = &DNTRare::pAAAR_;
    m_mapProbFunc["AARA_"] = &DNTRare::pAARA_;
    m_mapProbFunc["AARB_AB"] = &DNTRare::pAARB_AB;
    m_mapProbFunc["AARB_BA"] = &DNTRare::pAARB_BA;
    m_mapProbFunc["AABR_AB"] = &DNTRare::pAABR_AB;
    m_mapProbFunc["AABR_BA"] = &DNTRare::pAABR_BA;
    m_mapProbFunc["AARR"] = &DNTRare::pAARR;
    m_mapProbFunc["BARA"] = &DNTRare::pBARA;
    m_mapProbFunc["ABRA"] = &DNTRare::pABRA;
    m_mapProbFunc["BAAR"] = &DNTRare::pBAAR;
    m_mapProbFunc["ABAR"] = &DNTRare::pABAR;
    m_mapProbFunc["ABBR"] = &DNTRare::pABBR;
    m_mapProbFunc["BABR"] = &DNTRare::pBABR;
    m_mapProbFunc["ABRB"] = &DNTRare::pABRB;
    m_mapProbFunc["BARB"] = &DNTRare::pBARB;
    m_mapProbFunc["ABRC_ABC"] = &DNTRare::pABRC_ABC;
    m_mapProbFunc["ABRC_ACB"] = &DNTRare::pABRC_ACB;
    m_mapProbFunc["ABRC_CAB"] = &DNTRare::pABRC_CAB;
    m_mapProbFunc["ABCR_ABC"] = &DNTRare::pABCR_ABC;
    m_mapProbFunc["ABCR_ACB"] = &DNTRare::pABCR_ACB;
    m_mapProbFunc["ABCR_CAB"] = &DNTRare::pABCR_CAB;
    m_mapProbFunc["ABRR"] = &DNTRare::pABRR;
    m_mapProbFunc["ARAR"] = &DNTRare::pARAR;
    m_mapProbFunc["ARRA"] = &DNTRare::pARRA;
    m_mapProbFunc["RARA"] = &DNTRare::pRARA;
    m_mapProbFunc["ARBR_AB"] = &DNTRare::pARBR_AB;
    m_mapProbFunc["ARBR_BA"] = &DNTRare::pARBR_BA;
    m_mapProbFunc["ARRB_AB"] = &DNTRare::pARRB_AB;
    m_mapProbFunc["ARRB_BA"] = &DNTRare::pARRB_BA;
    m_mapProbFunc["RARB"] = &DNTRare::pRARB;
    m_mapProbFunc["ARRR"] = &DNTRare::pARRR;
    m_mapProbFunc["RARR"] = &DNTRare::pRARR;
    m_mapProbFunc["RRRR"] = &DNTRare::pRRRR;
    m_mapProbFunc["AAAA"] = &DNTRare::pAAAA;
    m_mapProbFunc["AAAB"] = &DNTRare::pAAAB;
    m_mapProbFunc["AABC"] = &DNTRare::pAABC;
    m_mapProbFunc["AABB"] = &DNTRare::pAABB;
    m_mapProbFunc["ABAB"] = &DNTRare::pABAB;
    m_mapProbFunc["ABAC"] = &DNTRare::pABAC;
    m_mapProbFunc["ABCD"] = &DNTRare::pABCD;
  }
  
  public:
    NumericVector prob(vector<string> vstrComb){
      
      int nCombs = vstrComb.size();
      NumericVector vResult;
      
      if(nCombs == 1 && vstrComb[0].compare("all") == 0){
        map<string, probFunc>::iterator i = m_mapProbFunc.begin();
        
        while(i != m_mapProbFunc.end()){
          probFunc f =  i->second; //m_mapProbFunc[i->first];
          vResult[i->first]  = (this->*f)();
          i++;
        }
      }else{
        for(int i = 0; i < nCombs; i++){
          map<string, probFunc>::iterator it = m_mapProbFunc.find(vstrComb[i]);
          
          if(it != m_mapProbFunc.end()){
            probFunc f =  m_mapProbFunc[vstrComb[i]];      
            vResult[vstrComb[i]]  = (this->*f)();
          }else{
            vResult[vstrComb[i]] = NA_REAL;
          }
        }
      }
      
      return vResult;
    }
};

#endif
