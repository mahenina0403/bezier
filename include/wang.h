// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"

using namespace std;

void convert_to_wang_ball(const vector<vec2> values, const vector<double> weights, vector<vec2> *WB, vector<double> *Ww, int n);
vec2 rationalWangBall(const vector<vec2> WB, const vector<double> Ww, double t);
double rcond_WangBall(const vector<vec2> WB, const vector<double> Ww, double t);
double rcond_WangBall2(const vector<vec2> WB, const vector<double> Ww, double t);