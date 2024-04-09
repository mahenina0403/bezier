#include "VS.h"

vec2 RationalVS(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s;
	double b = 1;
	vec2 N(0);
	double D;
	double tmp;

	if (t < 0.5){
		s  = t / (1.0 - t);
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

void gen_VS_data( vector<vec2> *values, vector<double> *weights, int n){
	int b = 1;
	double tmp;
	for (int i=1; i<n; i++){
		b = b*(n-i+1)/i;
		tmp = b*(*weights)[i];
		(*values)[i] *= tmp;
		(*weights)[i] = tmp;
	}
}

vec2 RationalVS2(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s;
	double b = 1;
	vec2 N(0);
	double D;
	double tmp;
	int ni;
	if (t < 0.5){
		s  = t / (1.0 - t);
		N = values[n];
		D = weights[n];
		for (int i=1; i<=n; i++){
			ni = n-i;
			N = N*s + values[ni];
			D = D*s + weights[ni];
		}
	}
	else{
		s  = (1 - t) / t;
		N = values[0];
		D = weights[0];
		for (int i=1; i<=n; i++){
			N = N*s + values[i];
			D = D*s + weights[i];
		}
	}
	
	return N/D;
}
