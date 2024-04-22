// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include "wang.h"

vector<vector<unsigned long long>> Uk(int n){
	double ceil_n2 = ceil(n/2);
	vector<vector<unsigned long long>> u(n+1);
	
	for (int k=0; k<=n; k++)
		u[k].resize(n+1);

	long long pow2 = 1;
	for (int i=0;i<=ceil_n2-1;i++){
		int Ni = n-2-2*i;
		for (int k=0; k<=n; k++){
			if (i > k)
				u[k][i] = 0;
			else if (n-2-2*i < 0)
				u[k][i] = 0;
			else if (k==i)
				u[k][i] = 1;
			else{
				u[k][i] = u[k-1][i] * (Ni+i-k+1);
				u[k][i] = u[k][i]/(k-i);
			}

		}
		for (int k=0; k<=n; k++){
			u[k][i] *= pow2;
			u[n-k][n-i] = u[k][i];
		}
		pow2 *= 2;
	}
	return u;
}

vec2 at_k(const vector<vec2> P, const vector<vector<unsigned long long>> u, int k){
	int n = P.size()-1;
	int i;

	vec2 A = vec2(0);
	
	for (i=0; i<=n; i++){
		A = A + u[k][i] * P[i];
	}
	return (A);
}

double at_k(const vector<double> P, const vector<vector<unsigned long long>> u, int k){
	int n = P.size()-1;
	int i;

	double A = 0;
	for (i=0; i<=n; i++){
		A = A + u[k][i] * P[i];
	}
	return (A);
}


void convert_to_wang_ball(const vector<vec2> values, const vector<double> weights, vector<vec2> *WB, vector<double> *Ww, int n){

	vector<vector<unsigned long long>> u = Uk(n);

	vector<vec2> Values(n+1);
	vector<double> Weights(n+1);
	double constant = 1;
    for (int i=0; i<=n; i++){
        Values[i] =  values[i] * weights[i];
        Weights[i] =  weights[i];

        (*WB)[i] = vec2(0);
        (*Ww)[i] =0;
    }

    double floor_n2 = floor(n/2);
	double ceil_n2 = ceil(n/2);

	unsigned long long binomial = 1;

	(*WB)[0] = Values[0];
	(*Ww)[0] = Weights[0];
	(*WB)[n] = Values[n];
	(*Ww)[n] = Weights[n];
	int k = 1;
	while(k <= n-k){
		binomial = binomial * (n-k+1);
		binomial = binomial / k;

		constant = constant / 2;
		(*WB)[k] = (binomial*Values[k]-at_k((*WB),u,k)) * constant;
		(*Ww)[k] = (binomial*Weights[k]-at_k((*Ww),u,k)) * constant;

		int K = n-k;

		if (K==k)
			break;

		(*WB)[K] = (binomial*Values[K]-at_k((*WB),u,K)) * constant;
		(*Ww)[K] = (binomial*Weights[K]-at_k((*Ww),u,K)) * constant;
		k++;
	}
	for (int i=0; i<=n; i++)
        (*WB)[i] = (*WB)[i]/(*Ww)[i];

}

vec2 rationalWangBall(const vector<vec2> WB, const vector<double> Ww, double t){

	int n = (WB).size()-1;

	if (n < 2){
		cout << "Error: Wang-Ball curve needs 3 points at least" << endl;
		exit(1);
	}
	vector<vec2> N(n+1);
	vector<double> D(n+1);
	for(int i=0;i<n+1;i++){
		N[i] = (WB)[i];
		D[i] = (Ww)[i];
	}

	double t1 = 1 - t;

	int k = n;
	int jump = 0;
	while (k>2){
		if (k%2==1){
			int n1 = (k-1)/2;
			int n2 = (k+1)/2;

			double a = t1*D[n1];
			double b = t*D[n2];
			if (jump == 0)
				jump = n2+1;
			
			D[n1] = a+b;
			N[n1] = (a*N[n1]+b*N[n2])/D[n1];
			
		}
		else{
			int n2 = k/2;
			if (jump==0)
				jump=n2+1;

			double a = t1*D[n2-1];
			double b = t*D[n2];
			double c = t1*D[n2];
			double d = t*D[jump];
	
			D[n2-1] = a+b;
			D[n2] = c+d;
	
			N[n2-1] = (a*N[n2-1]+b*N[n2])/D[n2-1];
			N[n2] = (c*N[n2]+d*N[jump])/D[n2];

			jump++;
			
		}
		k--;
	}

	double a = t1*D[0];
	double b = t*D[1];
	double c = t1*D[1];
	double d = t*D[n];

	double wq = a + b;
	double wr = c + d;

	double e = t1*wq;
	double f = t*wr;
	double ww = e + f;

	vec2 Q = (a*N[0] + b*N[1])/wq;
	vec2 V = (c*N[1] + d*N[n])/wr;
	vec2 W = (e*Q + f*V)/ww;

	return W;
}
