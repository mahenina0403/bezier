#include <vector>
#include <cmath>

#include "deCasteljau.h"

#define pi 3.14159265

vector<double> get_barycentric_values(vector<double> values, vector<double> weights, int n, double t);
vector<double> get_barycentric_weights(vector<double> weights, int n, double t);
double barycentric(vector<double> values, vector<double> weights, int n, double t);