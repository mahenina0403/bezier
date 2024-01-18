#include "lader.h"

mpreal lader(vector<mpreal> values, int n, mpreal t){

	mpreal s;
	mpreal b;
	mpreal N;
	mpreal sk;

	if (t < 1/2){
		s  = t/(1-t);
		sk = 1;
		b = 1;
		N = values[n];
		for (int i=1; i<=n; i++){
			sk = sk * s * (n-i+1)/i;
			// b = b * (n-i+1)/i;
			N = N + sk * values[n-i];
		}
	}
	else{
		s  = (1-t)/t;
		sk = 1;
		b = 1;
		N = values[0];
		for (int i=1; i<=n; i++){
			sk = sk * s * (n-i+1)/i;
			// b = b * (n-i+1)/i;
			N = N + sk * values[i];
		}
	}
	return N;
}

mpreal PolynomialLader(vector<mpreal> values, int n, mpreal t){

	mpreal s;
	mpreal b;
	mpreal N;

	if (t < 1/2){
		s  = pow(1-t,n);
	}
	else{
		s  = pow(t,n);
	}

	return s*lader(values, n, t);
}

mpreal RationalLader(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
	vector<mpreal> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	mpreal numerator = lader(R, n, t);
	mpreal denominator = lader(weights, n, t);

	return numerator/denominator;
}