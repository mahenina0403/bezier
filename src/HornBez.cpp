// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include "HornBez.h"

vec2 RationalHornBez(const vector<vec2> values, const vector<double> weights, int n, double t){
	double s = 1 - t;
	double tk = 1;
	double b = 1;
	vec2 N = s*values[0];
	double D = s * weights[0];

	for (int i=1; i<n; i++){
		tk = tk * t;
		N = (N+tk*values[i])*s;
		D = (D+tk*weights[i])*s;
	}

	tk = tk * t;
	N = N + tk*values[n];
	D = D + tk*weights[n];
	return N/D;
}