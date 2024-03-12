#ifndef LTBEZIER_H_INCLUDED
#define LTBEZIER_H_INCLUDED

#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <numbers>
#include <stdexcept>

using namespace std;

#define CHEBYSHEV 100
#define UNIFORM 101

template <class Type>
const Type Pi(){
    return atan((Type)1)*4;
}
// #define pi Pi<Type>()
//DECASTELJAU

template <class Type>
Type PolynomialDeCasteljau(vector<Type> values, int n, Type t){
    Type s  = 1 - t;

	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++)
			values[i] = s*values[i] + t*values[i+1];
	}

	return values[0];
}

template <class Type>
Type RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
    vector<Type> num(n+1);
	vector<Type> den(n+1);
	for(int i=0; i<=n; i++){
        num[i] = (Type) weights[i]*values[i];
        den[i] = (Type) weights[i];
	}
	Type s  = 1 - (Type) t;

    for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++){
			num[i] = (Type)s*num[i] + t*num[i+1];
			den[i] = (Type)s*den[i] + t*den[i+1];
		}
	}

	return num[0]/den[0];
}

template <class Type>
Type FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t){
	vector<Type> R(n+1);
	vector<Type> w(n+1);
	for(int i=0; i<=n; i++){
        R[i] = (Type) values[i];
        w[i] = (Type) weights[i];
	}

	Type s = 1 -(Type) t;

	Type c1, c2, tmp;
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


//HORNBEZ
template <class Type>
Type RationalHornBez(vector<double> values, vector<double> weights, int n, double t){
	Type s = 1 -(Type) t;
	Type tk = 1;
	Type b = 1;
	Type tmp = (Type) weights[0]*s;
	Type N = tmp*(Type)values[0]*s;
	Type D = tmp;

	for (int i=1; i<n; i++){
		tk *= t;
		b *= (n+1-i)/(Type)i;
		tmp = tk*b*(Type)weights[i];
		N = (N+tmp*(Type)values[i])*s;
		D = (D+tmp)*s;
	}
	// if (n > 0){
	tk *= t;
	tmp = tk*(Type)weights[n];
	N += tmp*(Type)values[n];
	D += tmp;
	// }

	return N/D;
}

//LINEAR GEOMETRIC
template <class Type>
Type linearGeometric(vector<double> values, vector<double> weights, int n, double t){
	Type h = 1;
	Type u = 1 - (Type)t;
	unsigned int n1 = n + 1;
	Type R = values[0];

	if(t <= 0.5) {
		u = t / u;
		for(int k=1; k<=n; k++){
			h *= u * (n1-k) * weights[k];
			h /= (k * (Type)weights[k-1] + h);
			R = (1 - (Type)h) * R + h * values[k];
		}
	}
	else {
		u /= t;
		for(int k = 1; k <= n; k++) {
			h *= (n1-k) * (Type)weights[k];
			h /= (k * u * weights[k-1] + h);
			R = (1 - (Type)h) * R + h * values[k];
		}
	}

	return R;
}

//LADER
template <class Type>
Type RationalLader(vector<double> values, vector<double> weights, int n, double t){
    Type s;
	Type b = 1;
	Type sk = 1;
    Type N;
    Type D;

	if (t < 0.5){
	    s = t / (1 - (Type)t);
		N = (Type) weights[0]*values[0];
		D = (Type) weights[0];
		for (int i=1; i<=n; i++){
			sk *= s * (n-i+1) / (Type)i;
			auto tmp = sk*weights[i];
			N += tmp*values[i];
			D += tmp;
		}
	}
	else{
	    s = (1 - (Type)t)/t;
		N = (Type) weights[n]*values[n];
		D = (Type) weights[n];
		for (int i=1; i<=n; i++){
			sk *= s * (n-i+1) / (Type)i;
			auto tmp = sk*weights[n-i];
			N += tmp*values[n-i];
			D += tmp;
		}
	}

	return N/D;
}

//VS
template <class Type>
Type RationalVS(vector<double> values, vector<double> weights, int n, double t){
	Type s;
	Type b = 1;
	Type N;
	Type D;

	if (t < 0.5){
		s  = t / (1 - (Type)t);
		N = (Type)weights[n]*values[n];
		D = (Type)weights[n];
		for (int i=1; i<=n; i++){
			b *= (n-i+1) / (Type)i;
			auto tmp = (Type)weights[n-i]*b;
			N = N*s + tmp*values[n-i];
			D = D*s + tmp;
		}
	}
	else{
		s  = (1 - (Type)t) / t;
		N = (Type)weights[0]*values[0];
		D = (Type)weights[0];
		for (int i=1; i<=n; i++){
			b *= (n-i+1) / (Type)i;
			auto tmp = (Type)weights[i]*b;
			N = N*s + tmp*values[i];
			D = D*s + tmp;
		}
	}

	return N/D;
}

//WANG
template <class Type>
Type n_choose_k(int n, int k){
	if (k==0)
		return 1;

	return n_choose_k<Type>(n,k-1) * (n-k+1) / (Type)k;
}

template <class Type>
Type polynomialWangBall(vector<Type> values, int n, double t){
	vector<Type> f(n);

	Type t1 = 1 - (Type)t;
	if (n >= 3){
		if (n%2 == 0){
			int n2 = n/2;
			for (int i=0; i<= (n2 - 2); i++)
				f[i] = values[i];
			f[n2-1] = values[n2-1]*t1 + (Type)t*values[n2];
			f[n2] = values[n2]*t1 + (Type)t*values[n2+1];
			for(int i=n2+1; i<=n-1; i++)
				f[i] = values[i+1];
		}else{
			for (int i=0; i<= (n-3)/2; i++)
				f[i] = values[i];
			f[(n-1)/2] = values[(n-1)/2]*t1 + (Type)t*values[(n+1)/2];
			for(int i=(n+1)/2; i<=n-1; i++)
				f[i] = values[i+1];
		}
		return polynomialWangBall<Type>(f, n-1, t);
	}
	f[0] = values[0]*t1 + (Type)t*values[1];
	f[1] = values[1]*t1 + (Type)t*values[2];
	f[0] = f[0]*t1 + (Type)t*f[1];

	return f[0];
}

template <class Type>
Type at_k(vector<Type> P, int k){
	int n = P.size()-1;
	int i;

	Type A = 0;
	Type B = 0;

	Type floor_n2 = floor((Type)n/2);
	Type ceil_n2 = ceil((Type)n/2);

	int start = n-k+1;
	int end = n-k;
	if (k==floor_n2)
		start = k+2;
	if (k==ceil_n2)
		end = k-2;

	if (k <= floor_n2){
		for (i=0; i<=k-1; i++)
			A = A + pow((Type)2,i) * n_choose_k<Type>(n-2-2*i, k-i) * P[i];
		for (i=start; i<=n;i++)
			B = B + pow((Type)2,n-i) * n_choose_k<Type>(2*i-2-n,i-k) * P[i];
	}
	// else if (k ==floor_n2){
	// 	for (i=0; i<k; i++)
	// 		A = A + pow((Type)2,i) * n_choose_k<Type>(n-2-2*i, k-i) * P[i];
	// 	for (i=k+2; i<n+1;i++)
	// 		B = B + pow((Type)2,n-i) * n_choose_k<Type>(2*i-2-n,i-k) * P[i];
	// }
	else if (k >= ceil_n2){
		for (i=0; i<=end; i++)
			A = A + pow((Type)2,i) * n_choose_k<Type>(n-2-2*i,k-i) * P[i];
		for (i=k+1; i<=n; i++)
			B = B + pow((Type)2,n-i) * n_choose_k<Type>(2*i-2-n,i-k) * P[i];
	}

	return (A+B)/n_choose_k<Type>(n,k);
}


template <class Type>
vector<vector<Type>> convert_to_wang_ball_stable(vector<double> values, vector<double> weights, int n){
	typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	vector<Type> Values(n+1);
	vector<Type> Weights(n+1);
    for (int i=0; i<=n; i++){
        Values[i] = (Type) values[i] * weights[i];
        Weights[i] = (Type) weights[i];
    }

    Type floor_n2 = floor((Type)n/2);
	Type ceil_n2 = ceil((Type)n/2);

	vector<vector<Type>> R(2);
	R[0].resize(n+1);
	R[1].resize(n+1);

	Type binomial = 1;
	R[0][0] = Values[0];
	R[1][0] = Weights[0];
	R[0][n] = Values[n];
	R[1][n] = Weights[n];
	int k = 1;
	while(k <= n-k){
		binomial = binomial * (n-k+1)/(Type)k;
		Type constant = 0;
		constant = binomial / pow((Type)2,k);
		R[0][k] = (Values[k]-at_k<Type>(R[0],k)) * constant;
		R[1][k] = (Weights[k]-at_k<Type>(R[1],k)) * constant;

		int K = n-k;

		if (K==k)
			break;

		constant = binomial / pow((Type)2,n-K);
		R[0][K] = (Values[K]-at_k<Type>(R[0],K)) * constant;
		R[1][K] = (Weights[K]-at_k<Type>(R[1],K)) * constant;

		// cout << k << ", " << K << endl;

		k++;
	}
	for (int i=0; i<=n; i++)
        R[0][i] = R[0][i]/R[1][i];

	return R;
}

template <class Type>
vector<vector<Type>> convert_to_wang_ball(vector<double> values, vector<double> weights, int n){
	typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	vector<vector<Type>> R(2);
	R[0].resize(n+1);
	R[1].resize(n+1);

	int i,j;

	Matrix data(2,n+1);
	for (i=0; i<= n; i++){
		data(0,i) = (Type)values[i]*weights[i];
		data(1,i) = (Type)weights[i];
	}

	Matrix M(n+1,n+1);
	if (n%2 == 0){
		auto n2 = n/2;
		for (i=0; i<=n; i++){
			for (j=0; j<=n; j++){
				if (i <= n2-1 && i <= j && j <= n-2-i){
					M(i,j) = pow((Type)2,i) * n_choose_k<Type>(n-2-2*i,j-i)/n_choose_k<Type>(n,j);
				}else if (n2+1 <= i && max(0, n-2-i) <= j && j <= i){
					M(i,j) = pow((Type)2,n-i) * n_choose_k<Type>(2*i-2-n, i-j)/n_choose_k<Type>(n,j);
				}else if (i==j && i==n2){
					M(i,j) = pow((Type)2,n2)/n_choose_k<Type>(n,n2);
				}else M(i,j) = 0;
			}
		}
	}else{
		for (i=0; i<=n; i++){
			for (j=0; j<=n; j++){
				if (i <= (n-3)/2 && i <= j && j <= n-2-i){
					M(i,j) = pow((Type)2,i) * n_choose_k<Type>(n-2-2*i,j-i)/n_choose_k<Type>(n,j);
				}else if ((n+3)/2 <= i && max(0, n-2-i) <= j && j <= i){
					M(i,j) = pow((Type)2,n-i) * n_choose_k<Type>(2*i-2-n, i-j)/n_choose_k<Type>(n,j);
				}else if (i==j && (i==(n-1)/2 || i==(n+1)/2)){
					M(i,j) = pow((Type)2,(n-1)/2)/n_choose_k<Type>(n,(n-1)/2);
				}else M(i,j) = 0;
			}
		}
	}

	M = M.inverse();
	data *= M;

	for (i=0; i<=n; i++){
		R[0][i] = data(0,i)/data(1,i);
		R[1][i] = data(1,i);
	}

	return R;
}

template <class Type>
Type rationalWangBall(vector<vector<Type>> R, int n, double t){

    Type numerator = polynomialWangBall<Type>(R[0], n, t);
	Type denominator = polynomialWangBall<Type>(R[1], n, t);

	return numerator/denominator;
}

// template <class Type>
// Type rationalWangBall_2(vector<vector<Type>> R, double t){

	// int n = R[0].size()-1;
	// vector<Type> N(n+1);
	// vector<Type> D(n+1);
	
	// N = R[0];
	// D = R[1];

	// Type t1 = 1 - (Type)t;
	// Type a,b,c,d,e,f,wq,wr,ww;
	// int k = n;
	// int k1,k2;
	// int shift;
	// int n_jump = 0;
	// while (k > 2){
	// 	if (k%2==1){
	// 		k1 = (int)(n-1)/2;
	// 		k2 = (int)(n+1)/2;
	// 		// if (n_jump==0)
	// 		// 	n_jump = k2;
	// 		a = t1*D[k1];
	// 		b = t*D[k2];
	// 		wq = a+b;
	// 		D[k1] = wq;
	// 		N[k1] = (a*N[k1]+b*N[k2])/wq;
	// 		shift = k2;
	// 		// n_jump++;
	// 		// cout << k1 << ", " << n_jump << endl;
	// 	}else{
	// 		k2 = (int)k/2;
	// 		// if (n_jump==0)
	// 		// 	n_jump = k2+1;
	// 		a = t1*D[k2-1];
	// 		b = t*D[k2];
	// 		wq = a+b;
	// 		D[k2-1] = wq;
	// 		N[k2-1] = (a*N[k2-1]+b*N[k2])/wq;
			
	// 		a = t1*D[k2];
	// 		b = t*D[n_jump];
	// 		wq = a+b;
	// 		D[k2] = wq;
	// 		N[k2] = (a*N[k2]+b*N[n_jump])/wq;
	// 		shift = k2+1;
	// 		// cout << k2 << ", " << n_jump << endl;
	// 	}
	// 	// for (int i=shift; i<k; i++){
	// 	// 	D[i] = D[i+1];
	// 	// 	N[i] = N[i+1];
	// 	// }
	// 	n_jump++;
	// 	k--;


	// }

	// 	a = t1*D[0];
	// 	b = t*D[1];
	// 	c = t1*D[1];
	// 	d = t*D[n];

	// 	wq = a + b;
	// 	wr = c + d;

	// 	e = t1*wq;
	// 	f = t*wr;
	// 	ww = e + f;

	// 	Type Q = (a*N[0] + b*N[1])/wq;
	// 	Type V = (c*N[1] + d*N[n])/wr;
	// 	Q = (e*Q + f*V)/ww;

	// return Q;
	
// }

template <class Type>
Type rationalWangBall_2(vector<vector<Type>> R, double t){

	int n = R[0].size()-1;

	if (n < 2){
		cout << "Error: Wang-Ball curve needs 3 points at least" << endl;
		exit(1);
	}
	vector<Type> N(n+1);
	vector<Type> D(n+1);
	for(int i=0;i<n+1;i++){
		N[i] = R[0][i];
		D[i] = R[1][i];
	}

	Type t1 = 1 - (Type)t;

	int k = n;
	int jump = 0;
	while (k>2){
		if (k%2==1){
			int n1 = (k-1)/2;
			int n2 = (k+1)/2;

			Type a = t1*D[n1];
			Type b = t*D[n2];
			if (jump == 0)
				jump = n2+1;
			
			D[n1] = a+b;
			N[n1] = (a*N[n1]+b*N[n2])/D[n1];
			
		}
		else{
			int n2 = k/2;
			if (jump==0)
				jump=n2+1;

			Type a = t1*D[n2-1];
			Type b = t*D[n2];
			Type c = t1*D[n2];
			Type d = t*D[jump];
	
			D[n2-1] = a+b;
			D[n2] = c+d;
	
			N[n2-1] = (a*N[n2-1]+b*N[n2])/D[n2-1];
			N[n2] = (c*N[n2]+d*N[jump])/D[n2];

			jump++;
			
		}
		k--;
	}

	Type a = t1*D[0];
	Type b = t*D[1];
	Type c = t1*D[1];
	Type d = t*D[n];

	Type wq = a + b;
	Type wr = c + d;

	Type e = t1*wq;
	Type f = t*wr;
	Type ww = e + f;

	Type Q = (a*N[0] + b*N[1])/wq;
	Type V = (c*N[1] + d*N[n])/wr;
	Type W = (e*Q + f*V)/ww;

	return W;
	
}


// template <class Type>
// Type rationalWangBall_2(vector<vector<Type>> R, double t){

// 	int n = R[0].size()-1;

// 	Type t1 = 1 - (Type)t;

// 	if (n==2){
// 		Type a = t1*R[1][0];
// 		Type b = t*R[1][1];
// 		Type c = t1*R[1][1];
// 		Type d = t*R[1][2];

// 		Type wq = a + b;
// 		Type wr = c + d;

// 		Type e = t1*wq;
// 		Type f = t*wr;
// 		Type ww = e + f;

// 		Type Q = (a*R[0][0] + b*R[0][1])/wq;
// 		Type V = (c*R[0][1] + d*R[0][2])/wr;
// 		Type W = (e*Q + f*V)/ww;

// 		return W;
// 	}

// 	vector<vector<Type>> S(2);
// 	S[0].resize(n);
// 	S[1].resize(n);

// 	if (n%2==1){
// 		int n1 = (n-1)/2;
// 		int n2 = (n+1)/2;

// 		for (int i=0; i<= (n-3)/2; i++){
// 			S[0][i] = R[0][i];
// 			S[1][i] = R[1][i];
// 		}
		
// 		Type a = t1*R[1][n1];
// 		Type b = t*R[1][n2];

// 		S[1][n1] = a+b;
// 		S[0][n1] = (a*R[0][n1]+b*R[0][n2])/S[1][n1];

// 		for (int i=n2; i<= n-1; i++){
// 			S[0][i] = R[0][i+1];
// 			S[1][i] = R[1][i+1];
// 		}
// 	}
// 	else{
// 		int n2 = n/2;

// 		for (int i=0; i < n/2; i++){
// 			S[0][i] = R[0][i];
// 			S[1][i] = R[1][i];
// 		}
		
// 		Type a = t1*R[1][n2-1];
// 		Type b = t*R[1][n2];
// 		Type c = t1*R[1][n2];
// 		Type d = t*R[1][n2+1];

// 		S[1][n2-1] = a+b;
// 		S[1][n2] = c+d;

// 		S[0][n2-1] = (a*R[0][n2-1]+b*R[0][n2])/S[1][n2-1];
// 		S[0][n2] = (c*R[0][n2]+d*R[0][n2+1])/S[1][n2];

// 		for (int i=n2+1; i<= n-1; i++){
// 			S[0][i] = R[0][i+1];
// 			S[1][i] = R[1][i+1];
// 		}
// 	}
//     return rationalWangBall_2<Type>(S, t);
	
// }

//BARYCENTRIC
template <class Type>
vector<vector<Type>> get_data(vector<double> values, vector<double> weights, int n, int distribution = UNIFORM){
	Type pi;
	if (distribution==CHEBYSHEV)
    	pi = Pi<Type>();


    vector<Type> Values(n+1);
    vector<Type> W(n+1);
    vector<Type> ones(n+1);
    for (int i=0; i<=n; i++){
        Values[i] = (Type) values[i] * weights[i];
        W[i] = (Type) weights[i];
        ones[i]=1;
    }

    Type t;
    vector<vector<Type>> N(2);
    N[0].resize(n+1);
    N[1].resize(n+1);

    Type binomial = 1;
    Type lagrange_weight;
    for (int i=0; i<=n; i++){
    	if (distribution==UNIFORM)
        	t = i/(Type)n;
        else if (distribution==CHEBYSHEV){
                t = (cos(i/(Type)n * pi)+1)/2;
                
                if (i==0 || i==n)
                    lagrange_weight = 0.5;
                else
                    lagrange_weight = 1;
        }
        Type zf = RationalVS<Type>(Values, ones, n, t);
        Type z = RationalVS<Type>(W, ones, n, t);
        
        if (distribution==UNIFORM){
        	N[0][i] = binomial * zf;
        	N[1][i] =  binomial * z;
        	binomial *= (n+1-(i+1))/(Type)(i+1);
        }
        else if (distribution == CHEBYSHEV){
            N[0][i] = lagrange_weight * zf;
            N[1][i] =  lagrange_weight * z;
        }
        
    }

    return N;
}

//BARYCENTRIC

// template <class Type>
// vector<Type> get_nodes(int n){
// 	const Type pi = Pi<Type>();
// 	vector<Type> T(n+1);
// 	for (int i=0; i<=n; i++)
// 		// T[i] = (cos(i/(Type)n * pi)+1)/2;
// 		T[i] = i/n;

// 	return T;
// }

// template <class Type>
// vector<Type> get_homogeneous_values(vector<double> values, vector<double> weights, vector<Type> T, int n){
//     // const Type pi = Pi<Type>();
//     vector<Type> Values(n+1);
//     for (int i=0; i<=n; i++)
//         Values[i] = (Type) values[i] * weights[i];

//     Type t;
//     vector<Type> N(n+1);
//     for (int i=0; i<=n; i++){
//         // t = cos(i/(Type)n * pi);
//         // t = (t+1)/2;
//         t = T[i];
//         N[i] = PolynomialDeCasteljau<Type>(Values, n, t);
//     }
//     return N;
// }

// template <class Type>
// vector<Type> get_barycentric_weights(vector<double> weights, vector<Type> T, int n){
//     // const Type pi = Pi<Type>();
//     vector<Type> W(n+1);
//     for(int i=0; i<=n; i++)
//         W[i] = (Type) weights[i];

//     vector<Type> N(n+1);
//     Type t;
//     Type lagrange_weight;
//     int sgn = 1;

//     for (int i=0; i<=n; i++){
//         // t = cos(i/(Type)n * pi);
//         // t = (t+1)/2;
//         t = T[i];
//         if (i==0 || i==n)
//             lagrange_weight = 0.5;
//         else
//             lagrange_weight = 1;

//         N[i] = sgn * lagrange_weight * PolynomialDeCasteljau<Type>(W, n, t);
//         sgn = -sgn;
//     }

//     return N;
// }

// template <class Type>
// void get_barycentric_data(vector<Type>& values, vector<Type> weights, int n){
//     int sgn = 1;
//     Type lagrange_weight;
//     for (int i=0; i<=n; i++){
//         if (i==0 || i==n)
//             lagrange_weight = 0.5;
//         else
//             lagrange_weight = 1;

//         values[i] *= sgn * lagrange_weight / weights[i];
//         sgn = -sgn;
//     }
// }

template <class Type>
Type barycentric(vector<Type> V, vector<Type> W, int n, double t, int distribution=UNIFORM){
	Type pi;
	if (distribution==CHEBYSHEV)
    	pi = Pi<Type>();

    //evaluation
    Type N = 0;
    Type D = 0;
    Type r;
    int sgn = 1;
    for (int i=0; i<=n; i++) {
        if (distribution==UNIFORM)
        	r = (Type)t-i/(Type)n;
        else if (distribution==CHEBYSHEV)
        	r = t - (cos(i/(Type)n * pi) + 1) / 2;	

        // Type r = (Type)t - T[i];
        if (r == 0)
            return V[i]/W[i];
        r = sgn / r;
        N += r * V[i];
        D += r * W[i];
        sgn = -sgn;
    }
    // for (int i=0)
    // cout << D << endl;
    return N/D;
}
#endif // LTBEZIER_H_INCLUDED