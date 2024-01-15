#include "VS.h"

double VS(vector<double> values, int n, double t){

	double s;
	double b;
	double N;

	if (t < 1/2){
		s  = t/(1-t);
		b = 1;
		N = values[n];
		for (int i=1; i<=n; i++){
			b = b * (n-i+1)/i;
			N = N*s + b * values[n-i];
		}
	}
	else{
		s  = (1-t)/t;
		b = 1;
		N = values[0];
		for (int i=1; i<=n; i++){
			b = b * (n-i+1)/i;
			N = N*s + b * values[i];
		}
	}

	return N;
}

double PolynomialVS(vector<double> values, int n, double t){

	double s;
	double b;
	double N;

	if (t < 1/2){
		s  = pow(1-t,n);
	}
	else{
		s  = pow(t,n);
	}

	return s*VS(values, n, t);
}

double RationalVS(vector<double> values, vector<double> weights, int n, double t){
	vector<double> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	double numerator = PolynomialVS(R, n, t);
	double denominator = PolynomialVS(weights, n, t);

	return numerator/denominator;
}