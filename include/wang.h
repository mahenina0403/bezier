#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Dense>

// #include <iostream>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

using namespace std;

vector<vector<double>> convert_to_wang_ball(vector<double> values, vector<double> weights, int n);
double rationalWangBall(vector<vector<double>> values, int n, double t);