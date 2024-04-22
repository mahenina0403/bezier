// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include <vector>
#include <cmath>
#include <iostream>
#include "vec2.h"
#include "VS.h"

using namespace std;

#define Pi() atan(1)*4

#define CHEBYSHEV 100
#define UNIFORM 101

vector<double> compute_nodes(int n,int distribution);
void get_data(const vector<vec2> values, const vector<double> weights, const vector<double> T, int n, vector<vec2> *Q, vector<double> *beta,int distrubution);
vec2 barycentric(const vector<vec2> V, const vector<double> W, const vector<double> T, int n, double t);
vector<vec2> barycentric_2(const vector<vec2> V, const vector<double> W, const vector<double> T, int n, double t);