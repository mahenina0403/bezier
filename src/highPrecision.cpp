// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include "highprecision.h"

mpreal RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
    vector<mpreal> num(n+1);
	vector<mpreal> den(n+1);
	for(int i=0; i<=n; i++){
        num[i] = (mpreal) weights[i]*values[i];
        den[i] = (mpreal) weights[i];
	}
	mpreal s  = 1 - (mpreal) t;

    for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++){
			num[i] = (mpreal)s*num[i] + t*num[i+1];
			den[i] = (mpreal)s*den[i] + t*den[i+1];
		}
	}

	return num[0]/den[0];
}

mpreal FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
	vector<mpreal> R(n+1);
	vector<mpreal> w(n+1);
	for(int i=0; i<=n; i++){
        R[i] = (mpreal) values[i];
        w[i] = (mpreal) weights[i];
	}

	mpreal s = 1 -(mpreal) t;

	mpreal c1, c2, tmp;
	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++) {
			c1 = s*w[i];
			c2 = (mpreal)t*w[i+1];
			tmp = c1 + c2;
			c1 = c1/tmp;
			c2 = 1 - c1;
			w[i] = tmp;
			R[i] = R[i]*c1 + c2*R[i+1];
		}
	}

	return R[0];
}