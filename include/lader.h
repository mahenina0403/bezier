#include <vector>
#include <cmath>
#include "mpreal.h"

using namespace std;
using namespace mpfr;

mpreal lader(vector<mpreal> values, int n, mpreal t);
mpreal PolynomialLader(vector<mpreal> values, int n, mpreal t);
mpreal RationalLader(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);