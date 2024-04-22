// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include <vector>
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

mpreal RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);
mpreal FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);