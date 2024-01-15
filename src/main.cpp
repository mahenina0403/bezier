#include "main.h"

int main() {
    int n = 3;
    vector<double> f;
    vector<double> beta;
    f.resize(n+1);
    beta.resize(n+1);

    f[0] = 1;
    f[1] = 1;
    f[2] = 1;
    f[3] = 1;

    beta[0] = 1;
    beta[1] = 1;
    beta[2] = 1;
    beta[3] = 1;

    double t = 0.5;
    cout << "Polynomial de Casteljau: " << PolynomialDeCasteljau(f,n,t) << endl;
    cout << "Rational de Casteljau: " << RationalDeCasteljau(f,beta,n,t) << endl;
    cout << "Rational (Farin) de Casteljau: " << FarinRationalDeCasteljau(f,beta,n,t) << endl;
    cout << "Polynomial VS: " << PolynomialVS(f,n,t) << endl;
    cout << "Rational VS: " << RationalVS(f,beta,n,t) << endl;
    cout << "Polynomial HornBez: " << PolynomialHornBez(f,n,t) << endl;
    cout << "Rational HornBez: " << RationalHornBez(f,beta,n,t) << endl;
    return 0;
}