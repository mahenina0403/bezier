#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include "mpreal.h"

using namespace std;
using namespace mpfr;

typedef Eigen::Matrix<mpreal, Eigen::Dynamic, Eigen::Dynamic> Matrix;

vector<vector<mpreal>> convert_to_wang_ball(vector<mpreal> values, vector<mpreal> weights, int n);
mpreal rationalWangBall(vector<vector<mpreal>> values, int n, mpreal t);