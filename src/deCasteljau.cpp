#include "deCasteljau.h"

mpreal PolynomialDeCasteljau(vector<mpreal> values, int n, mpreal t){
	mpreal s  = 1-t;

	vector<mpreal> R(n+1);
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

mpreal RationalDeCasteljau(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
	vector<mpreal> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	mpreal numerator = PolynomialDeCasteljau(R, n, t);
	mpreal denominator = PolynomialDeCasteljau(weights, n, t);

	return numerator/denominator;
}

mpreal FarinRationalDeCasteljau(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
	vector<mpreal> R(n+1);
	vector<mpreal> w(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = values[n-i];
		w[i] = weights[n-i];
	}

	mpreal s = 1-t;

	mpreal c1, c2, tmp;
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