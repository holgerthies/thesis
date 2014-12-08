/***
** Simulate Computations of Cray X-MP
** 
**/
#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include "dyadic_interval.h"
#include <vector>
#include <algorithm>
using namespace iRRAM;
using std::vector;


const DYADIC a = 3.6;
const DYADIC p0=0.4;

// computes the sqrt for an interval.
// Uses the real sqrt function to compute the end points
dyadic_interval sqrt(const dyadic_interval& x){
	dyadic_interval ans;
	DYADIC p1 = sqrt(REAL(x.left));
	DYADIC p2 = sqrt(REAL(x.right));
	ans.left = std::min(p1, p2);
	ans.right = std::max(p1,p2);
	return ans;
}

// computes the logistic map f(x) = a*x(1-x).
// Given the dyadic interval x outputs a dyadic_interval I  
// so that f(t) is in I for all t in x
dyadic_interval logistic_map(const dyadic_interval& x){
	dyadic_interval x_inv; // 1-x
	x_inv.left=1-x.right;
	x_inv.right = 1-x.left;
	return a*x*x_inv;
}

// computes one solution of the inverse logistic map of the dyadic_interval x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
dyadic_interval inverse_logistic_map(const dyadic_interval& x, bool left){
	if(left){
		//return 0.5-sqrt(0.25-x*1/a);
	}
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
	
	auto po = iRRAM_exec(pseudo_orbit, 1000, -10000);

	cout << setRwidth(100);
	cout << po[1000] << "\n";
	//for(DYADIC x : po)
	//	cout << x <<"\n";

}

