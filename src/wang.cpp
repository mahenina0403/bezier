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
	vector<vector<mpreal>> f;
	f.resize(n+1);
	for (int i=0; i<=n; i++){
		f[i].resize(n+1);
	}

	
	for (int i=0; i<=n; i++){
		f[0][i] = values[i];
	}

	mpreal t1 = 1-t;
	for (int k=n; k>=3; k--){
		if (k%2 == 0){
			int k2 = k/2;
			for (int i=0; i<= (k2 - 2); i++){
				f[n+1-k][i] = f[n-k][i];
			}
			f[n+1-k][k2-1] = f[n-k][k2-1]*t1 + t*f[n-k][k2];
			f[n+1-k][k2] = f[n-k][k2]*t1 + t*f[n-k][k2+1];
			for(int i=k2+1; i<=k-1; i++){
				f[n+1-k][i] = f[n-k][i+1];
			}
		}else{
			for (int i=0; i<= (k-3)/2; i++){
				f[n+1-k][i] = f[n-k][i];
			}
			f[n+1-k][(k-1)/2] = f[n-k][(k-1)/2]*t1 + t*f[n-k][(k+1)/2];
			for(int i=(k+1)/2; i<=k-1; i++){
				f[n+1-k][i] = f[n-k][i+1];
			}
		}
	}
	f[n-1][0] = f[n-2][0]*t1 + t*f[n-2][1];
	f[n-1][1] = f[n-2][1]*t1 + t*f[n-2][2];
	f[n][0] = f[n-1][0]*t1 + t*f[n-1][1];
	return f[n][0];
}

mpreal rationalWangBall(vector<vector<mpreal>> values, int n, mpreal t){
	mpreal numerator = polynomialWangBall(values[0], n, t);
	mpreal denominator = polynomialWangBall(values[1], n, t);

	// cout << numerator << "/" << denominator << "=";
	return numerator/denominator;
}