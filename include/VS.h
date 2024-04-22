// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"

using namespace std;

vec2 RationalVS(vector<vec2> values, vector<double> weights, int n, double t);
void gen_VS_data(const vector<vec2> f, const vector<double> beta, vector<vec2> *values, vector<double> *weights, int n);
vec2 RationalVS2(const vector<vec2> values, const vector<double> weights, int n, double t);
