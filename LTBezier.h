#ifndef LTBEZIER_H_INCLUDED
#define LTBEZIER_H_INCLUDED

#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <numbers>

using namespace std;

template <class Type>
const Type Pi(){
    return atan((Type)1)*4;
}

//DECASTELJAU

template <class Type>
Type PolynomialDeCasteljau(vector<Type> values, int n, Type t){
    Type s  = 1 - t;

	for (int j=1; j<=n; j++){
		for (int i=0; i<=n-j; i++)
			values[i] = t*values[i] + s*values[i+1];
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
			num[i] = (Type)t*num[i] + s*num[i+1];
			den[i] = (Type)t*den[i] + s*den[i+1];
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
			c1 = t*w[i];
			c2 = s*w[i+1];
			tmp = c1 + c2;
			c1 = c1/tmp;
			c2 = 1 - c1;
			w[i] = tmp;
			R[i] = R[i]*c1 + c2*R[i+1];
		}
	}

	return R[0];
}

//BARYCENTRIC
template <class Type>
vector<Type> get_homogeneous_values(vector<double> values, vector<double> weights, int n){
    const Type pi = Pi<Type>();
    vector<Type> Values(n+1);
    for (int i=0; i<=n; i++)
        Values[i] = (Type) values[i] * weights[i];

    Type t;
    vector<Type> N(n+1);
    for (int i=0; i<=n; i++){
        t = cos(i/(Type)n * pi);
        t = (t+1)/2;
        N[i] = PolynomialDeCasteljau<Type>(Values, n, t);
    }

    return N;
}

template <class Type>
vector<Type> get_barycentric_weights(vector<double> weights, int n){
    const Type pi = Pi<Type>();
    vector<Type> W(n+1);
    for(int i=0; i<=n; i++)
        W[i] = (Type) weights[i];

    vector<Type> N(n+1);
    Type t;
    Type lagrange_weight;
    int sgn = 1;

    for (int i=0; i<=n; i++){
        t = cos(i/(Type)n * pi);
        t = (t+1)/2;
        if (i==0 || i==n)
            lagrange_weight = 0.5;
        else
            lagrange_weight = 1;

        N[i] = sgn * lagrange_weight * PolynomialDeCasteljau<Type>(W, n, t);
        sgn = -sgn;
    }

    return N;
}

template <class Type>
void get_barycentric_data(vector<Type>& values, vector<Type> weights, int n){
    int sgn = 1;
    Type lagrange_weight;
    for (int i=0; i<=n; i++){
        if (i==0 || i==n)
            lagrange_weight = 0.5;
        else
            lagrange_weight = 1;

        values[i] *= sgn * lagrange_weight / weights[i];
        sgn = -sgn;
    }
}

template <class Type>
Type barycentric(vector<Type> V, vector<Type> W, int n, double t){
    const Type pi = Pi<Type>();

    //evaluation
    Type N = 0;
    Type D = 0;
    for (int i=0; i<=n; i++) {
        Type r = (Type)t - (cos(i/(Type)n * pi) + 1) / 2;
        if (r == 0)
            return V[i];
        r = W[i]/r;
        N += r * V[i];
        D += r;
    }

    return N/D;
}

//HORNBEZ
template <class Type>
Type RationalHornBez(vector<double> values, vector<double> weights, int n, double t){
	Type s = 1 -(Type) t;
	Type tk = 1;
	Type b = 1;
	Type N = (Type) weights[n]*values[n]*s;
	Type D = (Type) weights[n]*s;

	for (int i=1; i<n; i++){
		tk *= t;
		b *= (n+1-i)/(Type)i;
		N = (N+tk*b*weights[n-i]*values[n-i])*s;
		D = (D+tk*b*weights[n-i])*s;
	}
	if (n > 0){
		tk *= t;
		N += tk*weights[0]*values[0];
		D += tk*weights[0];
	}

	return N/D;
}

//LADER
template <class Type>
Type RationalLader(vector<double> values, vector<double> weights, int n, double t){
    Type s;
	Type b = 1;
	Type sk = 1;
    Type N;
    Type D;

	if (t < 1/2){
	    s = t / (1 - (Type)t);
		N = (Type) weights[n]*values[n];
		D = (Type) weights[n];
		for (int i=1; i<=n; i++){
			sk *= s * (n-i+1) / (Type)i;
			N += sk*weights[n-i]*values[n-i];
			D += sk*weights[n-i];
		}
	}
	else{
	    s = (1 - (Type)t)/t;
		N = (Type) weights[0]*values[0];
		D = (Type) weights[0];
		for (int i=1; i<=n; i++){
			sk *= s * (n-i+1) / (Type)i;
			N += sk*weights[i]*values[i];
			D += sk*weights[i];
		}
	}

	return N/D;
}

//LINEAR GEOMETRIC
template <class Type>
Type linearGeometric(vector<double> values, vector<double> weights, int n, double t){
	Type h = 1;
	Type u = 1 - (Type)t;
	unsigned int n1 = n + 1;
	Type R = values[n];

	if(t <= 0.5) {
		u = t / u;
		for(int k=1; k<=n; k++){
			h *= u * (n1-k) * weights[n-k];
			h /= (k * (Type)weights[n-k+1] + h);
			R = (1 - (Type)h) * R + h * values[n-k];
		}
	}
	else {
		u /= t;
		for(int k = 1; k <= n; k++) {
			h *= (n1-k) * (Type)weights[n-k];
			h /= (k * u * weights[n-k+1] + h);
			R = (1 - (Type)h) * R + h * values[n-k];
		}
	}

	return R;
}

//VS
template <class Type>
Type RationalVS(vector<double> values, vector<double> weights, int n, double t){
	Type s;
	Type b = 1;
	Type N;
	Type D;

	if (t < 1/2){
		s  = t / (1 - (Type)t);
		N = (Type)weights[0]*values[0];
		D = (Type)weights[0];
		for (int i=1; i<=n; i++){
			b *= (n-i+1) / (Type)i;
			N = N*s + (Type)weights[i]*values[i]*b;
			D = D*s + (Type)weights[i]*b;
		}
	}
	else{
		s  = (1 - (Type)t) / t;
		N = (Type)weights[n]*values[n];
		D = (Type)weights[n];
		for (int i=1; i<=n; i++){
			b *= (n-i+1) / (Type)i;
			N = N*s + (Type)weights[n-i]*values[n-i]*b;
			D = D*s + (Type)weights[n-i]*b;
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
vector<vector<Type>> convert_to_wang_ball(vector<double> values, vector<double> weights, int n){
	typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	vector<vector<Type>> R(2);
	R[0].resize(n+1);
	R[1].resize(n+1);

	int i,j;

	Matrix data(2,n+1);
	for (i=0; i<= n; i++){
		data(0,i) = (Type)values[n-i]*weights[n-i];
		data(1,i) = (Type)weights[n-i];
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
		R[0][i] = data(0,i);
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

#endif // LTBEZIER_H_INCLUDED
