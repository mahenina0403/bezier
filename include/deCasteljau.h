#include <vector>

#include "mpreal.h"

using namespace std;
using namespace mpfr;

mpreal PolynomialDeCasteljau(vector<mpreal> values, int n, mpreal t);
mpreal RationalDeCasteljau(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);
mpreal FarinRationalDeCasteljau(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);