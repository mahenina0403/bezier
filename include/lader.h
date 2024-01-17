#include <vector>
#include <cmath>
#include <mpfr.h>

using namespace std;

double lader(vector<double> values, int n, double t);
double PolynomialLader(vector<double> values, int n, double t);
double RationalLader(vector<double> values, vector<double> weights, int n, double t);