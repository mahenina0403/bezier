// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include <unsupported/Eigen/FFT>
#include <vector>
#include <complex>
#include <iostream>

#include "vec2.h"

using namespace std;

void toHomogeneousVector(const vector<vec2> values, const vector<double> weights, vector<double> *X, vector<double> *Y, int n);
vector<complex<double>> toFourier(const vector<complex<double>> values, int n);
vector<complex<double>> roots_of_unity(int n);
vec2 BernsteinFourrier(vector<complex<double>> X, vector<complex<double>> Y, vector<complex<double>> W, vector<complex<double>> r, int n, double t);
vector<vec2> BernsteinFourrier_2(const vector<complex<double>> X, const vector<complex<double>> Y, const vector<complex<double>> W, const vector<complex<double>> r, int n, double t);