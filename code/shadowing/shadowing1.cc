/***
** Simulate Computations of Cray X-MP
** 
**/
#include "iRRAM/lib.h"
#include "iRRAM/core.h"

using namespace iRRAM;

struct dyadic_interval{
	DYADIC left, right;
};


// computes the sum of two intervals
// [a,b]+[c,d] = [a+c, b+d]
DYADIC sum(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	ans.left = x.left+y.left;
	ans.right = x.right+y.right;
	return ans;
}

// computes the difference between two intervals
// [a,b]-[c,d] = [a-d, b-c]
DYADIC diff(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	ans.left = x.left-y.right;
	ans.right = x.right-y.left;
	return ans;
}

// computes the product between two intervals
// [a,b]*[c,d] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]
DYADIC mul(const dyadic_interval& x, const dyadic_interval& y){
	dyadic_interval ans;
	// pairwise product between all interval endpoints
	vector<DYADIC> pw_products = {x.left*y.left, x.left*y.right, x.right*y.left, y.right*y.right}; 
	ans.left = min_element(pw_products.begin(), pw_products.end());
	ans.right = max_element(pw_products.begin(), pw_products.end());
	return ans;
}

// Computes an upper bound for the distance of the 
// true orbit x_n to the noisy orbit p_n for N iterates
DYADIC shadowing_bound(DYADIC p0, int p){
	return approx(p0,p);
}


// computes the logistic map f(x) = a*x(1-x).
// Given the dyadic interval x outputs a dyadic_interval I  
// so that f(t) is in I for all t in x
dyadic_interval logistic_map(const dyadic_interval& x){

}


template DYADIC iRRAM_exec <DYADIC,int,DYADIC> 
(DYADIC (*) (DYADIC,int),DYADIC,int);
// main routine that internally calls iRRAM three times:
int main (int argc,char **argv)
{
	iRRAM_initialize(argc,argv);
	DYADIC p0=0.4;
	DYADIC d2;
	DYADIC a = 3.6;
	int p=-1000;

	d2=iRRAM_exec(shadowing_bound,p0,p);

	cout << setRwidth(100);
	cout << iRRAM_exec(shadowing_bound,p0,p)<<"\n";

}

