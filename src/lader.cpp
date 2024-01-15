#include "lader.h"

double lader(vector<double> values, int n, double t){

	double s;
	double b;
	double N;
	double sk;

	if (t < 1/2){
		s  = t/(1-t);
		sk = 1;
		b = 1;
		N = values[n];
		for (int i=1; i<=n; i++){
			sk = sk * s;
			b = b * (n-i+1)/i;
			N = N + sk * b * values[n-i];
		}
	}
	else{
		s  = (1-t)/t;
		sk = 1;
		b = 1;
		N = values[0];
		for (int i=1; i<=n; i++){
			sk = sk * s;
			b = b * (n-i+1)/i;
			N = N + sk * b * values[i];
		}
	}
	return N;
}

double PolynomialLader(vector<double> values, int n, double t){

	double s;
	double b;
	double N;

	if (t < 1/2){
		s  = pow(1-t,n);
	}
	else{
		s  = pow(t,n);
	}

	return s*lader(values, n, t);
}

double RationalLader(vector<double> values, vector<double> weights, int n, double t){
	vector<double> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	double numerator = PolynomialLader(R, n, t);
	double denominator = PolynomialLader(weights, n, t);

	return numerator/denominator;
}