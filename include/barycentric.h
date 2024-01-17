#include <vector>
#include <cmath>
#include "mpreal.h"
#include "deCasteljau.h"

#define Pi 3.14159265
using namespace mpfr;

vector<mpreal> get_homogeneous_values(vector<mpreal> values, vector<mpreal> weights, int n);
vector<mpreal> get_barycentric_weights(vector<mpreal> weights, int n);
vector<mpreal> get_barycentric_data(vector<mpreal> values, vector<mpreal> weights, int n);
mpreal barycentric(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t);