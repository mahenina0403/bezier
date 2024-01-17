#include <vector>
#include <cmath>
#include "mpreal.h"

using namespace std;
using namespace mpfr;

mpreal VS(vector<mpreal> values, int n, mpreal t);
mpreal PolynomialVS(vector<mpreal> values, int n, mpreal t);
mpreal RationalVS(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);