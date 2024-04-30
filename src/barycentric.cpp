// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include "barycentric.h"

double power(double x, int n){
	if (n==0)
		return 1;

	double x2 = power(x,int(n/2));
	if (n%2==0)
		return x2 * x2;
	else
		return x*x2*x2;
}

vector<double> compute_nodes(int n, int distribution){
	vector<double> T(n+1);

	double pi = Pi();

	for (int i=0; i<=n; i++){
		if (distribution==CHEBYSHEV)
			T[i] = (cos((double)i/n * pi)+1)/2;
		if (distribution==UNIFORM)
			T[i] = (double)i/n;
	}
	return T;
}

vector<double> get_at_ti(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s;
	vec2 N(0);
	double D;
	double tmp;

	double fact;
	double t1;
	if (t < 0.5){
		t1 = (1-t);
		s  = t / t1;
		fact = power(t1,n);
		N = values[n];
		D = weights[n];
		for (int i=1; i<=n; i++){
			int ni = n-i;
			N = N*s + values[ni];
			D = D*s + weights[ni];
		}
	}
	else{
		s  = (1 - t) / t;
		fact = power(t,n);
		N = weights[0]*values[0];
		D = weights[0];
		for (int i=1; i<=n; i++){
			N = N*s + values[i];
			D = D*s + weights[i];
		}
	}

	N = N/D;
	return {N.x(), N.y(), fact*D};

	// return {fact*N.x(), fact*N.y(), fact*D};
}

void get_data(const vector<vec2> values, const vector<double> weights, const vector<double> T, int n, vector<vec2> *Q, vector<double> *beta, int distribution){
	
	double sgn = 1;
	double t;
	
	double coeff = 1;
	vector<double> tmp;

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	gen_VS_data(values, weights, &g, &alpha, n);

	for (int i=0; i<=n; i++){
		t = T[i];
		tmp = get_at_ti(g, alpha, n, t);
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


vector<vec2> barycentric_2(const vector<vec2> V, const vector<double> W, const vector<double> T, int n, double t){
//evaluation
	vec2 N1 = vec2(0);
	vec2 N2 = vec2(0);
	double d1 = 0;
	double d2 = 0;
	double r;
	double u;

	int ni;
	for (int i=0; i<=n; i++) {
		ni = n-i;
		r = t-T[i];
		if (r==0)
			return {V[i], V[ni]};

		u = W[ni]/r;
		r = W[i]/r;
		N1 = N1 + r * V[i];
		d1 = d1 + r;

		N2 = N2 + u * V[ni];
		d2 = d2 + u;
	}
	return {N1/d1, N2/d2};
}

vec2 rcond_barycentric(const vector<vec2> V, const vector<double> W, const vector<double> T, int n, double t){
//evaluation
	vec2 N = vec2(0);
	vec2 D = vec2(0);
	double r;
	vec2 R = vec2(0);

	for (int i=0; i<=n; i++) {
		r = t-T[i];
		if (abs(r)<2e-3)
			return 1;
		R = W[i]*V[i]/r;
		N = N + R.abs();
		D = D + R;
	}

	return N%D.abs();
}