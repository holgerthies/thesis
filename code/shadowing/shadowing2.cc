/***
** Simulate Computations of Cray X-MP
** 
**/
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/multiprecision/mpfr.hpp>  // Defines the Backend type that wraps MPFR
#include "fixed_precision_interval.h"

namespace mp = boost::multiprecision;    

typedef mp::number<mp::mpfr_float_backend<30>>  clay_float;
typedef mp::number<mp::mpfr_float_backend<10>>  clay_double;

using std::cout;
using std::vector; 


/*
// thicken the given interval
// i.e. enlarge both sides of the interval by the 
// dyadic number given in the second parameter
dyadic_interval thicken(const dyadic_interval& x, const DYADIC& d){
	return dyadic_interval(x.left-d, x.right+d);
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
		return 0.5-sqrt(0.25-x*(1/a));
	}
	else{
		return 0.5+sqrt(0.25-x*(1/a));
	}
}*/


//  computes the logistic map for a dyadic number x
clay_float logistic_map(const clay_float& x, const clay_float& a){
	return a*x*(1-x);
}

/*
// compute the previous interval to a given interval
// using the algorithm from the paper
dyadic_interval prev_interval(dyadic_interval x, bool left){
	DYADIC thickening_constant = approx(REAL(pow(REAL(2),-90)),-90);
	//cout << thickening_constant << std::endl;
	DYADIC error_constant =approx(10*REAL(pow(REAL(2), -25)), -90);
	dyadic_interval x_thick = thicken(x, thickening_constant);
	cout << setRwidth(100);
	//cout << x_thick.left << " " << x_thick.right << std::endl;
	dyadic_interval x_prev = inverse_logistic_map(x, left);
	while(!logistic_map(x_prev).contains(x_thick)){
		x_prev = thicken(x_prev, thickening_constant);
	}
	//x_prev = thicken(x_prev, error_constant);
	return x_prev;
}*/

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
vector<clay_float> pseudo_orbit(const clay_float& p0, const clay_float& a, const int N){
	vector<clay_float> ans(N+1);
	ans[0] = p0;
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(ans[i-1],a);
	}
	return ans;
}
/*
// Computes an upper bound for the distance of the 
// shadow orbit x_n to the noisy orbit p_n for N iterates
DYADIC shadowing_bound(int N, int prec){
	//const int N=1000;
	vector<DYADIC> orbit = pseudo_orbit(N, prec);
	dyadic_interval x = dyadic_interval(orbit[N], orbit[N]);
	DYADIC beta = 0; // max distance between the p_n and the shadow orbit
	for(int i=N-1; i>=0; i--){
		bool left= (orbit[i] <= 0.5);
		x = prev_interval(x, left);
		beta = std::max(beta, x.max_dist(orbit[i]));
	}
	return beta;
}
*/


// main routine that internally calls iRRAM:
int main (int argc,char **argv)
{
	vector<clay_float> po = pseudo_orbit(0.4,3.6, 100);
	for(clay_float x : po)
		std::cout << x <<"\n";

}

