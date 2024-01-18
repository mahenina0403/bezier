#include "HornBez.h"

mpreal PolynomialHornBez(vector<mpreal> values, int n, mpreal t){

	mpreal s;
	mpreal tk;
	mpreal b;
	mpreal N;

	s = 1-t;
	tk = 1;
	b = 1;
	N = values[n]*s;

	for (int i=1; i<n; i++){
		tk = tk*t;
		b = b * (n+1-i)/i;
		N = (N+tk*b*values[n-i])*s;
	}
	if (n > 0){
		tk = tk*t;
		N = N + tk*values[0];
	}

	return N;
}

mpreal RationalHornBez(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
	vector<mpreal> R(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = weights[i]*values[i];
	}

	mpreal numerator = PolynomialHornBez(R, n, t);
	mpreal denominator = PolynomialHornBez(weights, n, t);

	return numerator/denominator;
}