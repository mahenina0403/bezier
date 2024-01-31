//gcc main.cpp -lstdc++ -lmpfr -o main.exe

#include <iostream>
#include <fstream>
#include "LTBezier.h"
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

int main(){
    int my_mpreal_precision = 1024;
    mpreal::set_default_prec(my_mpreal_precision);

    cout.precision(200);

    double t = 0.99;
    int n = 3;

	vector<double> f(n+1);
	vector<double> beta(n+1);

	f[0] = 0.1;
	f[1] = 0.00000002;
	f[2] = 0.0000001;
	f[3] = 0.3;

	beta[0] = 1;
	beta[1] = 0.1;
	beta[2] = 0.3;
	beta[3] = 1;

    mpreal deCast = RationalDeCasteljau<mpreal>(f,beta,n,t);
    mpreal farin = FarinRationalDeCasteljau<mpreal>(f,beta,n,t);
    mpreal vs = RationalVS<mpreal>(f,beta,n,t);
    mpreal hornBez = RationalHornBez<mpreal>(f,beta,n,t);
    mpreal lader = RationalLader<mpreal>(f,beta,n,t);
    mpreal linearGeo = linearGeometric<mpreal>(f,beta,n,t);

    vector<mpreal> V = get_homogeneous_values<mpreal>(f,beta,n);
    vector<mpreal> W = get_barycentric_weights<mpreal>(beta,n);
    get_barycentric_data<mpreal>(V,W,n);
    mpreal bary = barycentric<mpreal>(V,W,n,t);

    vector<vector<mpreal>> M = convert_to_wang_ball<mpreal>(f,beta,n);
    mpreal wangBall = rationalWangBall<mpreal>(M,n,t);

    cout << "DeCasteljau:      " << deCast << endl;
    cout << "FarinDeCasteljau: " << farin << endl;
    cout << "VS:               " << vs << endl;
    cout << "HornBez:          " << hornBez << endl;
    cout << "Lader:            " << lader << endl;
    cout << "linearGeometric:  " << linearGeo << endl;
    cout << "barycentric:      " << bary << endl;
    cout << "WangBall:         " << wangBall << endl;

}
