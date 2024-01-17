#include <vector>
#include <mpfr.h>

using namespace std;

double PolynomialDeCasteljau(vector<double> values, int n, double t);
double RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);
double FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);