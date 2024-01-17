#include <vector>
#include <mpfr.h>

using namespace std;

double PolynomialHornBez(vector<double> values, int n, double t);
double RationalHornBez(vector<double> values, vector<double> weights, int n, double t);