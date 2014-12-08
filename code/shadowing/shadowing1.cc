/***
** Simulate Computations of Cray X-MP
** 
**/
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <vector>
#include <algorithm>
using namespace iRRAM;
using std::vector;
using std::min_element;
using std::max_element;

const DYADIC a = 3.6;
const DYADIC p0=0.4;

struct dyadic_interval{
	DYADIC left, right;
};


// computes the sum of two intervals
// [a,b]+[c,d] = [a+c, b+d]
dyadic_interval sum(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	ans.left = x.left+y.left;
	ans.right = x.right+y.right;
	return ans;
}

// computes the difference between two intervals
// [a,b]-[c,d] = [a-d, b-c]
dyadic_interval diff(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	ans.left = x.left-y.right;
	ans.right = x.right-y.left;
	return ans;
}

// computes the product between two intervals
// [a,b]*[c,d] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]
dyadic_interval mul(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	// pairwise product between all interval endpoints
	vector<DYADIC> pw_products = {x.left*y.left, x.left*y.right, x.right*y.left, y.right*y.right}; 
	ans.left = *min_element(pw_products.begin(), pw_products.end());
	ans.right = *max_element(pw_products.begin(), pw_products.end());
	return ans;
}

// computes the scalar multiplication between a dyadic number d and an interval x
// d*[a, b] = [min(d*a, d*b), max(d*a, d*b)]
dyadic_interval scalar_mul(const DYADIC& d, const dyadic_interval& x){
	dyadic_interval ans;
	ans.left = std::min(d*x.left, d*x.right);
	ans.right = std::max(d*x.left, d*x.right);
	return ans;
}


// computes the logistic map f(x) = a*x(1-x).
// Given the dyadic interval x outputs a dyadic_interval I  
// so that f(t) is in I for all t in x
dyadic_interval logistic_map(const dyadic_interval& x){
	dyadic_interval x_inv; // 1-x
	x_inv.left=1-x.right;
	x_inv.right = 1-x.left;
	return scalar_mul(a, mul(x, x_inv));
}

//  computes the logistic map for a dyadic number x
REAL logistic_map(const REAL& x){
	return REAL(a)*x*(1-x);
}

// cuts all bits after 2^-prec
DYADIC force_approx(const REAL& x, const int prec){
	DYADIC ans = approx(x, prec);
	ans = INTEGER(scale(ans, -prec)); // multiply by 2^(prec) and round to integer
	ans = scale(ans, prec); // scale back
	return ans;
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The precision for all calculations is given by the input parameter p.
// The output is a vector of the result after each step
vector<DYADIC> pseudo_orbit(const int N, const int prec){
	vector<DYADIC> ans(N+1);
	ans[0] = p0;
	for(int i=1; i<=N; i++){
		ans[i] = force_approx(logistic_map(REAL(ans[i-1])),prec);
	}
	return ans;
}

// Computes an upper bound for the distance of the 
// true orbit x_n to the noisy orbit p_n for N iterates
DYADIC shadowing_bound(DYADIC p0, int p){
	dyadic_interval d;
	d.left = p0;
	d.right = p0;
	return approx(logistic_map(d).left,p);
}

template DYADIC iRRAM_exec <DYADIC,int,DYADIC> 
(DYADIC (*) (DYADIC,int),DYADIC,int);


template vector<DYADIC> iRRAM_exec <int,int,vector<DYADIC>> 
(vector<DYADIC> (*) (int,int),int,int);

// main routine that internally calls iRRAM:
int main (int argc,char **argv)
{
	iRRAM_initialize(argc,argv);
	DYADIC d2;
	
	int p=-1000;
	auto po = iRRAM_exec(pseudo_orbit, 1000, -10000);

	cout << setRwidth(100);
	cout << po[1000] << "\n";
	//for(DYADIC x : po)
	//	cout << x <<"\n";

}

