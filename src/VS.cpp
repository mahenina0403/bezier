#include "VS.h"

mpreal VS(vector<mpreal> values, int n, mpreal t){

	mpreal s;
	mpreal b;
	mpreal N;

	if (t < 1/2){
		s  = t/(1-t);
		b = 1;
		N = values[0];
		for (int i=1; i<=n; i++){
			b = b * (n-i+1)/i;
			N = N*s + b * values[i];
		}
	}
	else{
		s  = (1-t)/t;
		b = 1;
		N = values[n];
		for (int i=1; i<=n; i++){
			b = b * (n-i+1)/i;
			N = N*s + b * values[n-i];
		}
	}

	return N;
}

mpreal PolynomialVS(vector<mpreal> values, int n, mpreal t){

	mpreal s;
	mpreal b;
	mpreal N;

	if (t < 1/2){
		s  = pow(1-t,n);
	}
	else{
		s  = pow(t,n);
	}

	return s*VS(values, n, t);
}

mpreal RationalVS(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
	vector<mpreal> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	mpreal numerator = VS(R, n, t);
	mpreal denominator = VS(weights, n, t);

	return numerator/denominator;
}