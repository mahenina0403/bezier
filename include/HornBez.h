#include <vector>
#include "mpreal.h"

using namespace std;
using namespace mpfr;

mpreal PolynomialHornBez(vector<mpreal> values, int n, mpreal t);
mpreal RationalHornBez(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);