#include "BF.h"

complex<double> power(complex<double> x, int n){
	if (n==0)
		return 1;

	complex<double> x2 = power(x,int(n/2));
	if (n%2==0)
		return x2 * x2;
	else
		return x*x2*x2;
}

complex<double> power_arg(complex<double> x, int n){
	double theta = arg(x);
	double u = norm(x);

	return complex<double>(u*cos(n*theta),u*sin(n*theta));
}

void toHomogeneousVector(const vector<vec2> values, const vector<double> weights, vector<double> *X, vector<double> *Y, int n){
	for (int i=0; i<= n; i++){
		(*X)[i] = weights[i]*values[i].x();
		(*Y)[i] = weights[i]*values[i].y();
	}
}

vector<complex<double>> toFourier(const vector<complex<double>> values, int n){
	vector<complex<double>> U(n+1);

	Eigen::FFT<double> fft;
	fft.inv(U,values);

	return U;
}

vector<complex<double>> roots_of_unity(int n){
	vector<complex<double>> w(n+1);
	const double PI = std::acos(-1);

	for (int k = 0; k <= n; ++k) {
        double angle = -2 * PI * k / (n+1);
        w[k] = complex<double>(cos(angle), sin(angle));
    }

    return w;
}

double real_ab(complex<double> a, complex<double> b){
	return real(a)*real(b)-imag(b)*imag(a);
}

vec2 BernsteinFourrier(vector<complex<double>> X, vector<complex<double>> Y, vector<complex<double>> W, vector<complex<double>> r, int n, double t){
	complex<double> x, y, z;

	complex<double> c;
	complex<double> one(1,0);
	complex<double> two(2,0);
	complex<double> t_(t,0);

	int N;

	x = X[0];
	y = Y[0];
	z = W[0];

	if (n%2==0){
		N = n/2+1;
	}else{
		N = (n+1)/2;
		c = pow(1-2*t,n);
		x += c*X[N];
		y += c*Y[N];
		z += c*W[N];
	}
	for(int i=1; i<N; i++){
		c = one + t_*(r[i]-one);
		c = power(c,n);

		x = x + two*real(c*X[i]);
		y = y + two*real(c*Y[i]);
		z = z + two*real(c*W[i]);

	}

	return vec2(real(x/z), real(y/z));
}

vector<vec2> BernsteinFourrier_2(const vector<complex<double>> X, const vector<complex<double>> Y, const vector<complex<double>> W, const vector<complex<double>> r, int n, double t){
	// complex<double> xs, ys, zs;
	// complex<double> xt, yt, zt;
	
	double xs, ys, zs;
	double xt, yt, zt;
	

	complex<double> c;
	complex<double> one(1,0);
	complex<double> two(2,0);
	complex<double> t_(t,0);
	complex<double> t1(1-t,0);


	int N;

	xt = 0;
	yt = 0;
	zt = 0;

	xs = xt;
	ys = yt;
	zs = zt;

	if (n%2==0){
		N = n/2+1;
	}else{
		N = (n+1)/2;
	}
	for(int i=1; i<N; i++){
		c = r[i]*t_ + t1;
		c = power(c,n);
		// c = power_arg(c,n);
		// xt = xt + real(c*X[i]);
		// yt = yt + real(c*Y[i]);
		// zt = zt + real(c*W[i]);

		xt = xt + real_ab(c,X[i]);
		yt = yt + real_ab(c,Y[i]);
		zt = zt + real_ab(c,W[i]);

		c = conj(c)/r[i];
		// xs = xs + real(c*X[i]);
		// ys = ys + real(c*Y[i]);
		// zs = zs + real(c*W[i]);

		xs = xs + real_ab(c,X[i]);
		ys = ys + real_ab(c,Y[i]);
		zs = zs + real_ab(c,W[i]);

	}

	// xt = two*xt + X[0];
	// yt = two*yt + Y[0];
	// zt = two*zt + W[0];

	// xs = two*xs + X[0];
	// ys = two*ys + Y[0];
	// zs = two*zs + W[0];

	xt = 2*xt + real(X[0]);
	yt = 2*yt + real(Y[0]);
	zt = 2*zt + real(W[0]);

	xs = 2*xs + real(X[0]);
	ys = 2*ys + real(Y[0]);
	zs = 2*zs + real(W[0]);

	if (n%2==1){
		double u = pow(1-2*t,n);
		// c = power_arg(c,n);
		xt += u*real(X[N]);
		yt += u*real(Y[N]);
		zt += u*real(W[N]);

		xs -= u*real(X[N]);
		ys -= u*real(Y[N]);
		zs -= u*real(W[N]);

		// xt += real_ab(c,X[N]);
		// yt += real_ab(c,Y[N]);
		// zt += real_ab(c,W[N]);

		// xs -= real_ab(c,X[N]);
		// ys -= real_ab(c,Y[N]);
		// zs -= real_ab(c,W[N]);		
	}
	// return {vec2(real(xt/zt), real(yt/zt)), vec2(real(xs/zs), real(ys/zs))};
	return {vec2(xt/zt, yt/zt), vec2(xs/zs, ys/zs)};
}