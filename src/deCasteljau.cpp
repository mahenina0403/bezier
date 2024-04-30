// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include "deCasteljau.h"

vec2 RationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t){
    vector<vec2> num(n+1);
	vector<double> den(n+1);
	for(int i=0; i<=n; i++){
        num[i] = weights[i]*values[i];
        den[i] = weights[i];
	}
	double s  = 1 - t;

    for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++){
			num[i] = s*num[i] + t*num[i+1];
			den[i] = s*den[i] + t*den[i+1];
		}
	}

	return num[0]/den[0];
}

vec2 FarinRationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t){
	vector<vec2> R(n+1);
	vector<double> w(n+1);
	for (int i=0; i<=n; i++) {
		R[i] = values[i];
		w[i] = weights[i];
	}

	double s = 1-t;

	double c1, c2, tmp;
	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++) {
			c1 = s*w[i];
			c2 = t*w[i+1];
			tmp = c1 + c2;
			c1 = c1/tmp;
			c2 = 1 - c1;
			w[i] = tmp;
			R[i] = R[i]*c1 + c2*R[i+1];
		}
	}
	return R[0];
}

vec2 rcond_rdc(const vector<vec2> values, const vector<double> weights, int n, double t){
	vector<vec2> num_abs(n+1);
	vector<vec2> num(n+1);
	for(int i=0; i<=n; i++){
        num[i] = weights[i]*values[i];
        num_abs[i] = num[i].abs();
	}
	double s  = 1 - t;

    for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++){
			num[i] = s*num[i] + t*num[i+1];
			num_abs[i] = s*num_abs[i] + t*num_abs[i+1];
		}
	}

	return num_abs[0]%num[0].abs();
}
