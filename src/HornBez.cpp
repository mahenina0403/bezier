#include "HornBez.h"

double PolynomialHornBez(vector<double> values, int n, double t){

	double s;
	double tk;
	double b;
	double N;

	s = 1-t;
	tk = 1;
	b = 1;
	N = values[n]*s;

	for (int i=1; i<n; i++){
		tk = tk*t;
		b = b * (n+1-i)/i;
		N = (N+b*tk*values[n-i])*s;
	}
	if (n > 0){
		tk = tk*t;
		N = N + tk*values[0];
	}

	return N;
}

double RationalHornBez(vector<double> values, vector<double> weights, int n, double t){
	vector<double> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	double numerator = PolynomialHornBez(R, n, t);
	double denominator = PolynomialHornBez(weights, n, t);

	return numerator/denominator;
}