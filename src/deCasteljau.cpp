#include "deCasteljau.h"

double PolynomialDeCasteljau(vector<double> values, int n, double t){
	double s  = 1-t;

	vector<double> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = values[i];
	}
	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++) {
			R[i] = t*R[i] + s*R[i+1];
		}
	}
	return R[0];  
}

double RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
	vector<double> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	double numerator = PolynomialDeCasteljau(R, n, t);
	double denominator = PolynomialDeCasteljau(weights, n, t);

	return numerator/denominator;
}

double FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
	vector<double> R(n+1);
	vector<double> w(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = values[i];
		w[i] = weights[i];
	}

	double s = 1-t;

	double c1, c2, tmp;
	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++) {
			tmp = s*w[i] + t*w[i+1];
			c1 = s*w[i]/tmp;
			c2 = 1 - c1;
			w[i] = tmp;
			R[i] = R[i]*c1 + c2*R[i+1];
		}
	}
	return R[0];
}