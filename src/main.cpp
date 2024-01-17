#include "main.h"

int main(int argc, char *argv[]) {

	mpreal::set_default_prec(1024);
	

	int n = 3;
	vector<mpreal> f;
	vector<mpreal> beta;
	f.resize(n+1);
	beta.resize(n+1);

	f[0] = 0.1;
	f[1] = 0.00000002;
	f[2] = 0.0000001;
	f[3] = 0.1;

	beta[0] = 1;
	beta[1] = 0.1;
	beta[2] = 0.3;
	beta[3] = 1;

	mpreal t = 0.5;
	if (argc > 1) t = atof(argv[1]);

	auto rdc = RationalDeCasteljau(f,beta,n,t);
	// cout << "Precision: " << setprecision(50) << endl;
	cout << "rdc: " << rdc << endl;
	

	cout << "rdc - fdc: " << rdc - FarinRationalDeCasteljau(f,beta,n,t) << endl;
	cout << "rdc - rvs: " << rdc - RationalVS(f,beta,n,t) << endl;
	cout << "rdc - rhb: " << rdc - RationalHornBez(f,beta,n,t) << endl;
	cout << "rdc - rld: " << rdc - RationalLader(f,beta,n,t) << endl;
	cout << "rdc - ltg: " << rdc - linearGeometric(f,beta,n,t) << endl;

	auto bar_values = get_homogeneous_values(f,beta,n);
	auto bar_weights = get_barycentric_weights(beta,n);
	bar_values = get_barycentric_data(bar_values, bar_weights,n);
	cout << "rdc - bar: " << rdc - barycentric(bar_values,bar_weights,n,t) << endl;

	auto W = convert_to_wang_ball(f,beta,n);
	cout << "rdc - rwb: " << rdc - rationalWangBall(W,n,t) << endl;
	return 0;
}