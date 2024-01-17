#include <vector>
#include <cmath>
#include <mpfr.h>

using namespace std;

double VS(vector<double> values, int n, double t);
double PolynomialVS(vector<double> values, int n, double t);
double RationalVS(vector<double> values, vector<double> weights, int n, double t);