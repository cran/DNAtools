#include <Rcpp.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <sstream>

#include "multinomial.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace{
	// a little helper that should IMHO be standardized
	template<typename T>
	std::size_t make_hash(const T& v){
		return std::hash<T>()(v);
	}
	
	// adapted from boost::hash_combine
	void hash_combine(std::size_t& h, const std::size_t& v){
		h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
	}
	
	// hash any container
	template<typename T>
	struct hash_container{
		size_t operator()(const T& v) const{
			size_t h = 0;
			for( const auto& e : v ) {
				hash_combine(h, make_hash(e));
			}
			return h;
		}
	};
}

namespace std{
	// support for vector<T> if T is hashable
	// (the T... is a required trick if the vector has a non-standard allocator)
	template<>
	struct hash<IntegerVector> : hash_container<IntegerVector> {};
	
	struct equal_to_intvec : binary_function<IntegerVector, IntegerVector, bool> {
		bool operator() (const IntegerVector& x, const IntegerVector& y) const{
			if (x.size() != y.size()){
				return false;
			}
			
			int n = x.size();
			
			for (int i = 0; i < n; ++i) {
				if (x[i] != y[i]) {
					return false;
				}
			}
			
			return true;
		}
	};
	
	// the same for map<T,U> if T and U are hashable
	template<typename... T>
	struct hash<map<T...>> : hash_container<map<T...>> {};
	
	// simply add more containers as needed
}

class probsObj{
	class Alpha{
	  friend class probsObj;
		
		IntegerVector counts;
		map<int, int> groupSizes;
		int n;
		unsigned long long w;
  
  
		public:
		//default constructor
		Alpha(){
			n = 0;
		}
		
		Alpha(const IntegerVector& v){
			counts = v;
			setWeights();
			n = accumulate(counts.begin(), counts.end(), 0);
		}
		
		Alpha(const Alpha &a){
			counts = a.counts;
			groupSizes = a.groupSizes;
			n = a.n;
			w = a.w;
		}
		
		const Alpha& operator=(const Alpha& a){
			counts = a.counts;
			groupSizes = a.groupSizes;
			n = a.n;
			w = a.w;
			
			return *this;
		}
		
		int getDim(void){
			return (int)counts.size();
		}
		
		long getWeight(void){
			return w;
		}
		
		List toList(void){
			List alpha;
			
			alpha["a"] = NumericVector(counts.begin(), counts.end());
			alpha["w"] = w;
			
			return(alpha);
		}
		
		friend ostream &operator<<( ostream &output, const Alpha &a ){
			ostringstream oss;
			
			output << '(';
			
			IntegerVector::const_iterator it = a.counts.begin();
			
			while(it != a.counts.end()){
				output << *it;
				it++;
				output << (it == a.counts.end() ? ')' : ',');
			}
			
			output << ' ' << a.w;
			
			return output;
		}
		
		bool operator<(const Alpha& rhs){
			return this->n < rhs.n;
		}
		
		private:
		unsigned long long fact(unsigned long long n){
			switch(n){
				case 0:
				case 1:
					return 1;
				case 2:
					return 2;
				case 3:
					return 6;
				case 4:
					return 24;
				case 5:
					return 120;
				case 6:
					return 720;
				case 7:
					return 5040;
				case 8:
					return 40320;
				case 9:
					return 362880;
				case 10:
					return 3628800;
				default:
				  return n * fact(n - 1);
			}
		}
		  
	  unsigned long long multinomCoeff(IntegerVector x){
	    int nx = x.size();
	    multinomial::SVI v(nx);
	    int i;
	    
	    for(i = 0; i < nx; i++){
	      v.at(i) = x[i];
	    }
	    
	    // if(useDouble){
	    //   double u = multinomial::multi<double>(v);
	    //   NumericVector r = NumericVector::create(u);
	    //   
	    //   return r;
	    // }//else
	    unsigned long long u = multinomial::multi<unsigned long long>(v);
	    //NumericVector r = NumericVector::create(u);
	    
	    return u;
	  }
		
		void setWeights(void){
			IntegerVector::iterator it = counts.begin();
			n = 0;
			
			while(it != counts.end()){
				if(groupSizes.find(*it) == groupSizes.end()){
					groupSizes[*it] = 1;
				}else{
					groupSizes[*it]++;
				}
				n += *it;
				it++; 
			}
			
			//Rcpp::Rcout << counts << endl;
			w = multinomCoeff(counts);
			//Rcpp::Rcout <<  "w = " <<  w << endl;
			
			/*it = counts.begin();
			int p = 1;
			
			while(it != counts.end()){
				w /= fact(*it);
			  Rcpp::Rcout << "p" << p++ << " = " << (*it) << "!";
				it++;
			} */
			
			map<int, int>::iterator group = groupSizes.begin();
			//int d = 1;
			
			while(group != groupSizes.end()){
				if(group->second > 1){
					w /= fact(group->second);
				}
				group++;
			}
		//	Rcpp::Rcout << w << endl << endl;
		}
	}; // end class

private:
	NumericVector p;
	double m_dTheta;
	vector<Alpha> A;
	unordered_map<IntegerVector, double, std::hash<IntegerVector>, equal_to_intvec > lookup2;
	
private:
	struct row_greater
	{
		bool operator()( const IntegerVector& lx, const IntegerVector& rx ) const {
			int n = lx.size();
			
			for (int i = 0; i < n; ++i) {
				if (lx[i] < rx[i]) {
					return true;
				} else if (lx[i] > rx[i]) {
					return false;
				} else { // lx[i] == r[i]
					// continue next column
				}
			}
			
			return false; // strict weak ordering: http://www.sgi.com/tech/stl/StrictWeakOrdering.html, X < X: false
		}
	};
	
	IntegerMatrix updateAlpha_(IntegerVector& a){
		int na = a.size();
		IntegerMatrix aa(na-1,na-1);
		IntegerVector heada = head(a,na-1);
		for(int i=0; i<na-1; i++){
			aa(i,_) = heada;
			aa(i,i) += a[na-1];
		}
		return(aa);
	}
	
	bool rows_equal(IntegerVector& v1, IntegerVector& v2) {
		int n = v1.size();
		
		for (int i = 0; i < n; i++) {
			if (v1[i] != v2[i]) {
				return false;
			}
		}
  
		return true;
	}
	
	List matrix_table(IntegerMatrix& A) {
		int nr = A.nrow();
		//int nc = A.ncol(); not used
		
		std::vector<IntegerVector> rows(nr); // was nr
		
		for (int i = 0; i < nr; ++i) {
			IntegerVector v = A(i, _);
			std::sort(v.begin(), v.end());
			rows[i] = v;
		}
		
		std::sort(rows.begin(), rows.end(), row_greater());
		
		std::vector<IntegerVector> unique_rows;
		std::vector<int> unique_rows_count;
		int counter = 1;
		
		for (int i = 1; i < nr; ++i) {
			if (rows_equal(rows[i], rows[i - 1])) {
				counter++;
				continue;
			}
			
			unique_rows.push_back(rows[i - 1]);
			unique_rows_count.push_back(counter);
			counter = 1;
		}
		
		// Last rows
		unique_rows.push_back(rows[nr-1]);
		unique_rows_count.push_back(counter);
		
		List ret;
		ret["rows"] = wrap(unique_rows);
		ret["counts"] = wrap(unique_rows_count);
		
		return ret;
	}
	
	string toString(const IntegerVector& a){
	  ostringstream oss;
	  IntegerVector::const_iterator i = a.begin();
	  while(i != a.end()){
	    oss << *i << ',';
	    i++;
	  }
	  return oss.str();
	}
	
	string toString(const IntegerVector& a, const IntegerVector& b){
	  ostringstream oss;
	  oss << toString(a) << ';' << toString(b);
	  return oss.str();
	}
		
	double Sa_(IntegerVector a){
		auto iter = lookup2.find(a);
	  
		if(iter != lookup2.end()){
	    return iter->second;
		}
		
		if (a.size() == 1){
			return(sum(pow(p, a[0])));
		}
		else{
			IntegerMatrix update_a = updateAlpha_(a);
			List perm = matrix_table(update_a);
			//List perm = matrix_table_weight(update_a, rep(1, update_a.nrow()));

			double res = 0;
			List permrow = perm["rows"];
			IntegerVector permcount = perm["counts"];
			
			for(int i = 0; i < permcount.size(); ++i){
				res += permcount[i] * Sa_(permrow[i]);
			}
			
			double theResult = Sa_(tail(a, 1)) * Sa_(head(a, a.size() - 1)) - res;
			
			/*
			if(theResult < 0){
			  Rcpp::Rcout << "update_a:" << endl;
			  Rcpp::print(update_a);
			  Rcpp::Rcout << endl << "perm:" << endl;
			  Rcpp::print(perm);
			  Rcpp::Rcout << endl << "permrow:" << endl;
			  Rcpp::print(permrow);
			  Rcpp::Rcout << endl << "permcount:" << endl;
			  Rcpp::print(permcount);
			  Rcpp::Rcout << endl << "a:" << endl;
			  Rcpp::print(a);
			  Rcpp::Rcout << endl << "theResult:" << endl;
			  Rcpp::Rcout << theResult << endl;
			}
			*/
			
			lookup2[a] = theResult;
			return(theResult);
		}
	}
	
	IntegerVector ek_(IntegerVector& x, int k){
		IntegerVector y(x.size());
	
		y[k - 1] = 1;
		return y;
	}
	
	double Sab_(IntegerVector a, IntegerVector b){
		int lenB = b.size();
		
		if (abs(m_dTheta) < 1e-12) { // m_dTheta == 0
			return  Sa_(a + b);
		}
		
		if(lenB == 0){
			return Sa_(a);
		}
		else if(is_true(any(b == 0))){
			return  Sab_(a, b[b != 0]);
		}
		else{
			if(b[lenB - 1] == 1 && lenB == 1){
			  return Sab_(a + ek_(a, lenB), head(b, lenB - 1));
			}
			else if(b[lenB - 1] == 1){
				 return (1.0 - m_dTheta) * Sab_(a + ek_(a, lenB), head(b, lenB - 1));
			}
			else{
			  return (b[lenB - 1] - 1) * m_dTheta * Sab_(a, b - ek_(b, lenB)) +
						   (1.0 - m_dTheta) * Sab_(a + ek_(a, lenB), b - ek_(b, lenB));
			}
		}
	}
	
	double Sab__(IntegerVector& b){
		if (abs(m_dTheta) < 1e-12) { // m_dTheta == 0
			return Sa_(b);
		}else{
			IntegerVector a(b.size());
			return(Sab_(a, b));
		}
	}
	
public:
	probsObj(void){
		m_dTheta = 0;
	}
	
	probsObj(double theta){
		m_dTheta = theta;
	}
	
	probsObj(NumericVector locusFreqs, double theta){
		m_dTheta = theta;
		p = locusFreqs;
	}
  
  void generateCompositions(int n){
    // Adapted from Python code and algorithm by
    // Jerome Kelleher (c) 2009
    // http://jeromekelleher.net/partitions.php
    vector<int> a(n+1);
    
    fill(a.begin(), a.end(), 0);
    
    int k = 1;
    int y = n - 1;
    
    while(k != 0){
      int x = a[k - 1] + 1;
      k -= 1;
      while(2*x <= y){
        a[k] = x;
        y -= x;
        k += 1;
      }
      int l = k + 1;
      while(x <= y){
        a[k] = x;
        a[l] = y;
        /*for(int i = 0; i < k + 2; i++)
        cout << a[i] << ' '; */
        A.push_back(IntegerVector(a.begin(), a.begin() + (k + 2)));
        //cout << endl;
        x += 1;
        y -= 1;
      }
      a[k] = x + y;
      y = x + y - 1;
      /*for(int i = 0; i < k + 1; i++)
      cout << a[i] << ' ';
      cout << endl;*/
      A.push_back(IntegerVector(a.begin(), a.begin() + k + 1));
    }
  }
  
  void printFreqs(){
    int nAlleles = p.size();
    
    for(int i = 0; i < nAlleles; i++){
      Rprintf("%d\t%.7f\n", i + 1, p[i]);
    }
  }
  
  List getCompositions(int numContrib){
    A.clear();
    generateCompositions(2 * numContrib);
    
    List comps;
    auto a = A.begin();
    
    while(a != A.end()){
      comps.push_back(a->toList());
      a++;
    }
    
    return comps;
  }
  
	NumericVector calcProbs(int numContrib, bool bZeroNegs = true){
		A.clear();
		generateCompositions(2 * numContrib);
		
		vector<Alpha>::iterator it = A.begin();

		// for each vector compute the probability for the alpha vector/composition
		NumericVector sums(2 * numContrib); // I believe these are initialized to zero
		
		/* 
		 * NOTE: 
		 * Sab__() can become negative with less than 2*m alleles;
		 * this is handled later.
		 */
		
		while(it != A.end()){
			sums[(*it).getDim() - 1] += (*it).w * Sab__((*it).counts);
		  it++;
		}
		
		double denominator = 1;
		int nAlleles = p.size();
		
		for(int i = 1; i <= 2 * (numContrib - 1); i++){
		  denominator *= 1 + i * m_dTheta;
		}
		
		// divide the probabilities by the denominator
		if(bZeroNegs){
		  bool isNeg = false;
		  for(int i = 0; i < 2 * numContrib; i++){
		    if(!isNeg && sums[i] < 0)
		      isNeg = true;
		    
		    sums[i] = ((isNeg || i >= nAlleles) ? 0 : sums[i] / denominator);
		  }
		  return sums;
		}else{
		  return sums / denominator;
		}
	}
  
  NumericMatrix calcProbs(int numContrib, List loci, bool bZeroNegs = true){
    A.clear();
    generateCompositions(2 * numContrib);
    int nLoci = loci.size();
		
		NumericMatrix result(nLoci, 2 * numContrib);
		
		Rcpp::CharacterVector rn = loci.names();
    Rcpp::CharacterVector cn(2 * numContrib);
    for(int i = 0; i < 2 * numContrib; i++){
      cn(i) = std::to_string(i + 1);
    }
		result.attr("dimnames") = Rcpp::List::create(rn, cn);
		
		double denominator = 1;
		for(int i = 1; i <= 2 * (numContrib - 1); i++)
		  denominator *= 1 + i * m_dTheta;
		
		auto	locus = loci.begin();
		int loc = 0;
		
		while(locus != loci.end()){
			// set the locus probabilities
			p = *locus;
		  int nAlleles = p.size();

		  // clear the lookup table
		 	lookup2.clear();
		  
			
			// iterate over combinations
			
			auto it = A.begin();

			// for each vector compute the probability for the alpha vector/composition

			while(it != A.end()){
			  result(loc, it->getDim() - 1) += it->w * Sab__(it->counts);
				it++;
			}
			
			for(int i = nAlleles; i < 2 * numContrib; i++)
			  result(loc, i) = 0;
			
			// divide the probabilities by the denominator
			if(bZeroNegs){
  			bool isNeg = false;
  			for(int i = 0; i < 2 * numContrib; i++){
  			  if(!isNeg && result(loc, i) < 0)
  			    isNeg = true;
  			  
  			  result(loc, i) = ((isNeg || i >= nAlleles) ? 0 : result(loc, i) / denominator);
  			}
			}else{
			  for(int i = 0; i < 2 * numContrib; i++){
			    result(loc, i) /= denominator;
			  }
			}
			
			loc++;
			locus++;
		}
		
		return result;
  }
  
  void resetProbs(NumericVector &newProbs){
    p = newProbs;
  }
  
  string printLookup(){
    auto i = lookup2.begin();
    ostringstream oss;
    
    while(i != lookup2.end()){
      oss << toString((i->first)) << '\t' << i->second << endl;
      i++;
    }
    
    return oss.str();
  }
};

