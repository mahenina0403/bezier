#include "VS.h"

vec2 RationalVS(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s;
	double b = 1;
	vec2 N;
	double D;
	double tmp;

	if (t < 0.5){
		s  = t / (1 - t);
		N = weights[n]*values[n];
		D = weights[n];
		for (int i=1; i<=n; i++){
			b = b*(n-i+1) / i;
			tmp = weights[n-i]*b;
			N = N*s + tmp*values[n-i];
			D = D*s + tmp;
		}
	}
	else{
		s  = (1 - t) / t;
		N = weights[0]*values[0];
		D = weights[0];
		for (int i=1; i<=n; i++){
			b = b*(n-i+1) / i;
			tmp = weights[i]*b;
			N = N*s + tmp*values[i];
			D = D*s + tmp;

		}
	}
	
	return N/D;
}