#include <vector>
#include <iostream>
#include <cmath>
#include "vec2.h"
using namespace std;

vec2 RationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t);
vec2 FarinRationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t);