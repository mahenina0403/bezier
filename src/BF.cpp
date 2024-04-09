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

vec2 BernsteinFourrier(vector<complex<double>> X, vector<complex<double>> Y, vector<complex<double>> W, vector<complex<double>> r, int n, double t){
	complex<double> x;
	complex<double> y;
	complex<double> z;

	complex<double> c;
	complex<double> unit(1,0);
	complex<double> t_(t,0);

	x = 0;
	y = 0;
	z = 0;

	for(int i=0; i<=n; i++){
		// c = complex<double>(1,0)+complex<double>(t,0)*(r[i]-complex<double>(1,0));
		c = unit + t_*(r[i]-unit);
		c = pow(c,n);

		x = x + c*X[i];
		y = y + c*Y[i];
		z = z + c*W[i];

	}

	return vec2(real(x/z), real(y/z));
}

vector<vec2> BernsteinFourrier_2(vector<complex<double>> X, vector<complex<double>> Y, vector<complex<double>> W, vector<complex<double>> r, int n, double t){
	complex<double> xs;
	complex<double> ys;
	complex<double> zs;

	complex<double> xt;
	complex<double> yt;
	complex<double> zt;

	complex<double> c;
	xs = 0;
	ys = 0;
	zs = 0;

	xt = 0;
	yt = 0;
	zt = 0;

	double s = 1-t;
	for(int i=0; i<=n; i++){
		c = complex<double>(1,0)+complex<double>(t,0)*(r[i]-complex<double>(1,0));
		c = pow(c,n);

		xt = xt + c*X[i];
		yt = yt + c*Y[i];
		zt = zt + c*W[i];

		c = c / conj(r[i]);

		xs = xs + c*conj(X[i]);
		ys = ys + c*conj(Y[i]);
		zs = zs + c*conj(W[i]);		
	}

	return {vec2(real(xt/zt), real(yt/zt)), vec2(real(xs/zs), real(ys/zs))};
}