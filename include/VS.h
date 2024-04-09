#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"

using namespace std;

vec2 RationalVS(vector<vec2> values, vector<double> weights, int n, double t);
// vector<vec2> RationalVS2(const vector<vec2> values, const vector<double> weights, int n, double t);

void gen_VS_data( vector<vec2> *values, vector<double> *weights, int n);
vec2 RationalVS2(const vector<vec2> values, const vector<double> weights, int n, double t);
