// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"

using namespace std;

void convert_to_wang_ball(const vector<vec2> values, const vector<double> weights, vector<vec2> *WB, vector<double> *Ww, int n);
vec2 rationalWangBall(const vector<vec2> WB, const vector<double> Ww, double t);
// double rcond_WangBall(const vector<vec2> WB, const vector<double> Ww, double t);
vec2 rcond_WangBall(const vector<vec2> WB, const vector<double> Ww, double t);