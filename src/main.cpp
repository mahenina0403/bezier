//gcc main.cpp -lstdc++ -lmpfr -o main.exe

#include <iostream>
#include <fstream>
#include <ctime>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "LTBezier.h"
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

void relative_error(vector<double> f, vector<double> beta, int n, double t){
    mpreal absolute = RationalDeCasteljau<mpreal>(f,beta,n,t);
    double decast = RationalDeCasteljau<double>(f,beta,n,t);
    double farin = FarinRationalDeCasteljau<double>(f,beta,n,t);
    double vs = RationalVS<double>(f,beta,n,t);
    double hornBez = RationalHornBez<double>(f,beta,n,t);
    double lader = RationalLader<double>(f,beta,n,t);
    double linearGeo = linearGeometric<double>(f,beta,n,t);

    // vector<double> T = get_nodes<double>(n);
    // vector<double> V = get_homogeneous_values<double>(f,beta,T,n);
    // vector<double> W = get_barycentric_weights<double>(beta,T,n);
    // get_barycentric_data<double>(V,W,n);
    // double bary = barycentric<double>(V,W,T,n,t);

    vector<vector<double>> D = get_data<double>(f,beta,n);
    // for(int i=0; i<=n; i++){
    //     cout << D[0][i] << ", " << D[1][i] << endl;
    // }
    double bary = barycentric<double>(D[0],D[1],n,t);

    vector<vector<double>>  M = convert_to_wang_ball_stable<double>(f,beta,n);
    double wangBall = rationalWangBall_2<double>(M,t);

    cout.precision(5);
    cout << "DeCasteljau:      " << absolute - (mpreal)decast << endl;
    cout << "FarinDeCasteljau: " << absolute - (mpreal)farin << endl;
    cout << "VS:               " << absolute - (mpreal)vs << endl;
    cout << "HornBez:          " << absolute - (mpreal)hornBez << endl;
    cout << "Lader:            " << absolute - (mpreal)lader << endl;
    cout << "linearGeometric:  " << absolute - (mpreal)linearGeo << endl;
    cout << "barycentric:      " << absolute - (mpreal)bary << endl;
    cout << "WangBall:         " << absolute - (mpreal)wangBall << endl;
}

void efficiency_comparison(vector<double> f, vector<double> beta, int n, double t, int Sample=1000){
    // int Sample = 1000;
    cout.precision(10);
    
    clock_t startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        RationalDeCasteljau<double>(f,beta,n,t);
    }
    clock_t endTime = clock();
    double t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "DeCasteljau:      " << t1 << "." << endl;

    startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        FarinRationalDeCasteljau<double>(f,beta,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "FarinDeCasteljau: " << t1 << "." << endl;
    
    startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        RationalVS<double>(f,beta,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "VS:               " << t1 << "." << endl;

    startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        RationalHornBez<double>(f,beta,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "HornBez:          " << t1 << "." << endl;

    startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        RationalLader<double>(f,beta,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "Lader:            " << t1 << "." << endl;

    startTime = clock();
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        linearGeometric<double>(f,beta,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "linearGeometric:  " << t1 << "." << endl;

    startTime = clock();
    vector<vector<double>> M = convert_to_wang_ball_stable<double>(f,beta,n);
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        rationalWangBall_2<double>(M,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "WangBall:         " << t1 << "." << endl;

    
    // vector<double> T = get_nodes<double>(n);
    // vector<double> V = get_homogeneous_values<double>(f,beta,T,n);
    // vector<double> W = get_barycentric_weights<double>(beta,T,n);
    // get_barycentric_data<double>(V,W,n);
    vector<double> T(n+1);
    // for (int i=0; i<=n; i++)
    //     T[i] = i/n;
    
    // T = get_nodes<double>(n);

    startTime = clock();
    vector<vector<double>> D = get_data<double>(f,beta,n);
    vector<double> V = D[0];
    vector<double> W = D[1];
    
    for (int i=0; i<=Sample; i++){
        t = i / (double) Sample;
        barycentric<double>(V,W,n,t);
    }
    endTime = clock();
    t1 = (endTime-startTime) / (double) CLOCKS_PER_SEC;
    cout << "barycentric:      " << t1 << "." << endl;
}

int main(int argc, char* argv[]){
    int my_mpreal_precision = 1024;
    mpreal::set_default_prec(my_mpreal_precision);

    cout.precision(200);

    double t = 0.5;
    int n = 20;

	vector<double> f(n+1);
	vector<double> beta(n+1);

    for(int i=0; i<=n; i++){
        f[i] = i*pow(-1,i);
        beta[i] = 1;
    }
    beta[0] = 2;
    beta[n] = 2;

    if (argc == 2 && strcmp(argv[1],"-r")==0)
        relative_error(f,beta,n,t);
    else if (argc == 3 && strcmp(argv[1],"-r")==0){
        t = atof(argv[2]);
        relative_error(f,beta,n,t);
    }
    else if (argc == 2 && strcmp(argv[1],"-e")==0)
        efficiency_comparison(f,beta,n,t);
    else if (argc == 3 && strcmp(argv[1],"-e")==0){
        efficiency_comparison(f,beta,n,t,atoi(argv[2]));
    }
    else{
        cout << "Command:" << endl;
        cout << setw(25) << left << "./compare -e";
        cout << "compare the efficiency with 1000 sample" << endl;
        cout << setw(25) << "./commpare -e [sample]";
        cout << "compare the efficiency with a user defined number of sample" << endl;
        cout << setw(25) << "./compare -r";
        cout << "ompare the relative error at t = 0.5" << endl;
        cout << setw(25) << "./compare -r [t]";
        cout << "compare the relative error at a user defined parameter t" << endl;
    }

    return 1;
}
