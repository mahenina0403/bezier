#include "main.h"

int main(int argc, char *argv[]) {
    int n = 3;
    vector<double> f;
    vector<double> beta;
    f.resize(n+1);
    beta.resize(n+1);

    f[0] = 1;
    f[1] = 1;
    f[2] = 1;
    f[3] = 1.3;

    beta[0] = 1;
    beta[1] = 1;
    beta[2] = 1;
    beta[3] = 1;

    double t = 0.5;
    if (argc > 1) t = atof(argv[1]);

    cout << "Polynomial de Casteljau: " << PolynomialDeCasteljau(f,n,t) << endl;
    cout << "Rational de Casteljau: " << RationalDeCasteljau(f,beta,n,t) << endl;
    cout << "Rational (Farin) de Casteljau: " << FarinRationalDeCasteljau(f,beta,n,t) << endl;
    cout << "Polynomial VS: " << PolynomialVS(f,n,t) << endl;
    cout << "Rational VS: " << RationalVS(f,beta,n,t) << endl;
    cout << "Polynomial HornBez: " << PolynomialHornBez(f,n,t) << endl;
    cout << "Rational HornBez: " << RationalHornBez(f,beta,n,t) << endl;
    cout << "Polynomial lader: " << PolynomialLader(f,n,t) << endl;
    cout << "Rational lader: " << RationalLader(f,beta,n,t) << endl;
    cout << "Linear-time geometric: " << linearGeometric(f,beta,n,t) << endl;

    auto bar_values = get_barycentric_values(f,beta,n);
    auto bar_weights = get_barycentric_weights(beta,n);
    cout << "Barycentric: " << barycentric(bar_values,bar_weights,n,t) << endl;

    auto W = convert_to_wang_ball(f,beta,n);
    // for (auto data: W){
    //     for (auto item: data){
    //         cout << item << ", ";
    //     }
    //     cout << endl;
    // }
    cout << "Wang-Ball: " << rationalWangBall(W,n,t) << endl;
    return 0;
}