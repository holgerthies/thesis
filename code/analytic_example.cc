#include "iRRAM.h"
#include "Functions/Analytic/BA_ANA.h"
using namespace iRRAM;

REAL factorial(const int n){
	if (n==0)
		return REAL(1);
	return factorial(n-1)*REAL(n);
}

REAL sinseries(const int n){
	if (n%2 == 0) return 0;
	if (n % 4 == 3) return -1/factorial(n); 
	return 1/factorial(n); 
}

void compute() {
	POWERSERIES<REAL> ps_sin(sinseries); 
	BA_ANA<REAL> sinus(ps_sin, 2,2);
	sinus.differentiate(2);
	rwrite(sinus(0.2), 20);
	iRRAM::cout << endl;
}
