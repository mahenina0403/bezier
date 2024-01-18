#include "wang.h"


mpreal n_choose_k(int n, int k){
	if (k==0){
		return 1;
	}
	return n_choose_k(n,k-1)*(n-k+1)/k;
}

vector<vector<mpreal>> convert_to_wang_ball(vector<mpreal> values, vector<mpreal> weights, int n){
	vector<vector<mpreal>> R;
	R.resize(2);
	R[0].resize(n+1);
	R[1].resize(n+1);

	int i,j;

	Matrix data;
	data.resize(2,n+1);
	for (i=0; i<= n; i++){
		data(0,i) = values[n-i]*weights[n-i];
		data(1,i) = weights[n-i];
	}
	

	Matrix M;
	M.resize(n+1,n+1);

	
	if (n%2==0){
		auto n2 = n/2;
		for (i=0; i<=n; i++){
			for (j=0; j<=n; j++){
				if (i <= n2-1 && i <= j && j <= n-2-i){
					M(i,j) = pow(2,i) * n_choose_k(n-2-2*i,j-i)/n_choose_k(n,j);
				}else if (n2+1 <= i && max(0, n-2-i) <= j && j <= i){
					M(i,j) = pow(2,n-i) * n_choose_k(2*i-2-n, i-j)/n_choose_k(n,j);
				}else if (i==j && i==n2){
					M(i,j) = pow(2,n2)/n_choose_k(n,n2);
				}else M(i,j) = 0;
			}
		}
	}else{
		for (i=0; i<=n; i++){
			for (j=0; j<=n; j++){
				if (i <= (n-3)/2 && i <= j && j <= n-2-i){
					M(i,j) = pow(2,i) * n_choose_k(n-2-2*i,j-i)/n_choose_k(n,j);
				}else if ((n+3)/2 <= i && max(0, n-2-i) <= j && j <= i){
					M(i,j) = pow(2,n-i) * n_choose_k(2*i-2-n, i-j)/n_choose_k(n,j);
				}else if (i==j && (i==(n-1)/2 || i==(n+1)/2)){
					M(i,j) = pow(2,(n-1)/2)/n_choose_k(n,(n-1)/2);
				}else M(i,j) = 0;
			}
		}
	}

	M = M.inverse();
	data = data*M;

	for (i=0; i<=n; i++){
		R[0][i] = data(0,i);
		R[1][i] = data(1,i);
	}
	return R;
}

mpreal polynomialWangBall(vector<mpreal> values, int n, mpreal t){
	vector<mpreal> f;
	f.resize(n);
	
	mpreal t1 = 1-t;
	if (n >= 3){
		if (n%2 == 0){
			int n2 = n/2;
			for (int i=0; i<= (n2 - 2); i++){
				f[i] = values[i];
			}
			f[n2-1] = values[n2-1]*t1 + t*values[n2];
			f[n2] = values[n2]*t1 + t*values[n2+1];
			for(int i=n2+1; i<=n-1; i++){
				f[i] = values[i+1];
			}
		}else{
			for (int i=0; i<= (n-3)/2; i++){
				f[i] = values[i];
			}
			f[(n-1)/2] = values[(n-1)/2]*t1 + t*values[(n+1)/2];
			for(int i=(n+1)/2; i<=n-1; i++){
				f[i] = values[i+1];
			}
		}
		return polynomialWangBall(f, n-1, t);
	}
	f[0] = values[0]*t1 + t*values[1];
	f[1] = values[1]*t1 + t*values[2];
	f[0] = f[0]*t1 + t*f[1];
	return f[0];
}

mpreal rationalWangBall(vector<vector<mpreal>> values, int n, mpreal t){
	mpreal numerator = polynomialWangBall(values[0], n, t);
	mpreal denominator = polynomialWangBall(values[1], n, t);

	// cout << numerator << "/" << denominator << "=";
	return numerator/denominator;
}