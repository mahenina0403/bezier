#include "main.h"

int main(int argc, char *argv[]) {

	int n = 3;
	vector<double> f;
	vector<double> beta;
	f.resize(n+1);
	beta.resize(n+1);

	f[0] = 1;
	f[1] = 0.2;
	f[2] = 0.1;
	f[3] = 1.1;

	beta[0] = 1;
	beta[1] = 0.1;
	beta[2] = 0.3;
	beta[3] = 1;

	double t = 0.5;
	if (argc > 1) t = atof(argv[1]);

	cout << "Rational de Casteljau: " << RationalDeCasteljau(f,beta,n,t) << endl;
	cout << "Rational (Farin) de Casteljau: " << FarinRationalDeCasteljau(f,beta,n,t) << endl;
	cout << "Rational VS: " << RationalVS(f,beta,n,t) << endl;
	cout << "Rational HornBez: " << RationalHornBez(f,beta,n,t) << endl;
	cout << "Rational lader: " << RationalLader(f,beta,n,t) << endl;
	cout << "Linear-time geometric: " << linearGeometric(f,beta,n,t) << endl;

	auto bar_values = get_homogeneous_values(f,beta,n);
	auto bar_weights = get_barycentric_weights(beta,n);
	bar_values = get_barycentric_data(bar_values, bar_weights,n);
	cout << "Barycentric: " << barycentric(bar_values,bar_weights,n,t) << endl;

	auto W = convert_to_wang_ball(f,beta,n);
	cout << "Wang-Ball: " << rationalWangBall(W,n,t) << endl;
	return 0;
}