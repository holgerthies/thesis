#include "iRRAM.h"
#include <vector>
#include <algorithm>


const float a=3.8;
const float p0=0.4;
using namespace iRRAM;
using std::endl;
using std::vector;
using std::max;

//  computes the logistic map for a dyadic number x
float logistic_map(const float x){
	return a*x*(1-x);
}

//  computes the logistic map for a dyadic number x
REAL logistic_map(const REAL& x){
	return REAL(a)*x*(REAL(1)-x);
}

// computes one solution of the inverse logistic map of the REAL x
// Which solution is chosen is decided by the parameter left. 
// If left is through the solution left of 0.5 is chosen.
REAL inverse_logistic_map(const REAL& x, bool left){
	if(left){
		return REAL(0.5)-sqrt(REAL(0.25)-x*1/REAL(a));
	}
	else{
		return REAL(0.5)+sqrt(REAL(0.25)-x*1/REAL(a));
	}
}

// computes a pseudo orbit (noisy orbit).
// Applies the logistic map N times. 
// The output is a vector of the result after each step
vector<float> pseudo_orbit(const int N){
	vector<float> ans(N+1);
	ans[0] = p0;
	for(int i=1; i<=N; i++){
		ans[i] = logistic_map(ans[i-1]);
	}
	return ans;
} 


void compute(){
	const int N=1000000;
	// max distance between pseudo orbit and closest of the real orbits currently computed
	REAL max_dist = 1000.0;
	vector<float> orbit = pseudo_orbit(N);
	REAL grid_size = power(REAL(2), -45);
	REAL left = REAL(orbit[N])-grid_size;
	REAL best_p = orbit[0];
	int n=2;
	while(!bound(max_dist, -18)){
		cout << n << endl;
		for(int i=0; i<n; i++){
			// max distance between real and pseudo orbit for current start point
			REAL start_point_dist=0.0; 
			// start point
			REAL p = left+i*grid_size;
			for(int j=N-1;j>=0; j--){
				bool left= (orbit[j] <= 0.5);
				p = inverse_logistic_map(p, left);
				REAL dist = abs(p-REAL(orbit[j]));
				// multivalued version of dist > start_point_dist
				if(positive(dist-start_point_dist,-100))
					start_point_dist = dist;
				if(!bound(dist, -15)) 
					break;
			}
			// multivalued version of max_dist > start_point_dist
			if(positive(max_dist-start_point_dist,-100)){
				best_p = p; // this is now p0
				max_dist = start_point_dist;
			}
		}
		cout<<best_p<<" "<<orbit[0]<<endl;
		grid_size /= 2;
		n *= 2;
		rshow(max_dist, 500);
		cout << endl;
	}
	cout<<best_p<<endl;
	REAL x=best_p;
	REAL md = 0.0;
	for(int i=0;i<=N;i++){
		REAL dist = abs(x-REAL(orbit[i]));
		if(positive(dist-md,-100))
			md = dist;
		x = logistic_map(x);
	}
	cout<<md<<endl;
	/*
	best_p=orbit[0];
	grid_size = power(REAL(2), -20);
	left = REAL(orbit[0])-grid_size;
	n=3;
	while(!bound(max_dist, -11)){
		cout << n << endl;
		for(int i=1; i<n; i++){
			REAL start_point_dist = 0.0;
			REAL p0=left+i*grid_size;
			REAL p=p0;
			for(int j=0;j<=N;j++){
				REAL dist = abs(p-REAL(orbit[j]));
				if(positive(dist-start_point_dist,-100))
					start_point_dist = dist;
			//cout << p << " "<<orbit[j]<<endl;
				p = logistic_map(p);
				if(positive(dist-power(REAL(2),-10),-100)){
					//cout<<"upto"<<j<<endl;
					break;
				}

			}
			if(positive(max_dist-start_point_dist,-100)){
				best_p = p0; // this is now p0
				max_dist = start_point_dist;
			}
			//return 0;
		}
		cout<<best_p<<" "<<orbit[0]<<endl;
		grid_size /= 2;
		n *= 2;
		rshow(max_dist, 500);
		cout << endl;
	}*/
}