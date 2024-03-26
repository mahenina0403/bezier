#include "HornBez.h"

vec2 RationalHornBez(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s = 1 - t;
	double tk = 1;
	double b = 1;
	double tmp = s * weights[0];
	vec2 N = tmp*values[0];
	double D = tmp;

	for (int i=1; i<n; i++){
		tk = tk * t;
		b = b*(n+1-i)/i;
		tmp = tk*b*weights[i];
		N = (N+tmp*values[i])*s;
		D = (D+tmp)*s;
	}

	tk = tk * t;
	tmp = tk*weights[n];
	N = N + tmp*values[n];
	D = D + tmp;
	return N/D;
}