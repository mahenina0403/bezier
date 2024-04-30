// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) [2024] [Chiara Fuda, Andriamahenina Ramanantoanina]


#include <vector>
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

mpreal RationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);
mpreal FarinRationalDeCasteljau(vector<double> values, vector<double> weights, int n, double t);