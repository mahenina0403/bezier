#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"

using namespace std;




double n_choose_k(int n, int k);
vector<vector<double>> uk(int n);
vec2 at_k(const vector<vec2> P, const vector<vector<double>> u, int k);
void convert_to_wang_ball(const vector<vec2> values, const vector<double> weights, vector<vec2> *WB, vector<double> *Ww, int n);
vec2 rationalWangBall(const vector<vec2> WB, const vector<double> Ww, double t);