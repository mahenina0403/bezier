#include "barycentric.h"

vector<double> compute_nodes(int n, int distribution){
	vector<double> T(n+1);

	double pi = Pi();

	for (int i=0; i<=n; i++){
		if (distribution==CHEBYSHEV)
			T[i] = (cos(i/n * pi)+1)/2;
		if (distribution==UNIFORM)
			T[i] = i/n;
	}
	return T;
}

vector<double> get_at_ti(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s;
	double b = 1;
	vec2 N;
	double D;
	double tmp;

	double fact;
	double t1;
	if (t < 0.5){
		s  = t / (1 - t);
		fact = 1;
		t1 = (1-t);
		N = weights[n]*values[n];
		D = weights[n];
		for (int i=1; i<=n; i++){
			b = b*(n-i+1) / i;
			tmp = weights[n-i]*b;
			N = N*s + tmp*values[n-i];
			D = D*s + tmp;
			fact = fact * t1;
		}
	}
	else{
		s  = (1 - t) / t;
		fact = 1;
		t1 = t;
		N = weights[0]*values[0];
		D = weights[0];
		for (int i=1; i<=n; i++){
			b = b*(n-i+1) / i;
			tmp = weights[i]*b;
			N = N*s + tmp*values[i];
			D = D*s + tmp;
			fact = fact * t1;
		}
	}
	N = N/D;
	return {N.x(), N.y(), fact*D};
}

void get_data(const vector<vec2> values, const vector<double> weights, const vector<double> T, int n, vector<vec2> *Q, vector<double> *beta, int distribution){
	
	double sgn = 1;
	double t;
	
	double coeff = 1;
	vector<double> tmp;

	for (int i=0; i<=n; i++){
		t = T[i];
		tmp = get_at_ti(values, weights, n, t);
		(*Q)[i] = vec2(tmp[0],tmp[1]);

		if (distribution==CHEBYSHEV){
			coeff = 1;
			if (i==0 || i==n)
        		coeff = 0.5;
		}
		(*beta)[i] =  sgn*coeff * tmp[2];

		if (distribution==UNIFORM)
			coeff = coeff * (n+1-(i+1))/(i+1);
		sgn = -sgn;
	}
}

vec2 barycentric(const vector<vec2> V, const vector<double> W, const vector<double> T, int n, double t){
//evaluation
	vec2 N = vec2(0);
	double D = 0;
	double r;
	int sgn = 1;
	for (int i=0; i<=n; i++) {
		r = t-T[i];
		if (r==0)
			return V[i];
		r = W[i] / r;
		N = N + r * V[i];
		D = D + r;
	}
	return N/D;
}