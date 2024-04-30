// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include <vector>
#include <iostream>
#include <cmath>
#include "vec2.h"
using namespace std;

vec2 RationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t);
vec2 FarinRationalDeCasteljau(const vector<vec2> values, const vector<double> weights, int n, double t);
vec2 rcond_rdc(const vector<vec2> values, const vector<double> weights, int n, double t);