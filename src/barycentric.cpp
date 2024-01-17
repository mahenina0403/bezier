#include "barycentric.h"

vector<double> get_barycentric_values(vector<double> values, vector<double> weights, int n){
  vector<double> Values;
  double t;
  Values.resize(n+1);

  for (int i=0; i<=n; i++){
    Values[i] = values[i] * weights[i];
  }

  vector<double> N;
  N.resize(n+1);
  for (int i=0; i<=n; i++){
    t = cos(Pi*i/n);
    N[i] = PolynomialDeCasteljau(Values, n, t);
  }
  return N;
}

vector<double> get_barycentric_weights(vector<double> weights, int n){
  vector<double> N;
  double t;
  double lagrange_weight;
  int sgn = -1;
  
  N.resize(n+1);
  for (int i=0; i<=n; i++){
    t = cos(Pi*i/n);
    if (i==0 || i==n){
      lagrange_weight = 0.5;
    }else{
      lagrange_weight = 1;
    }
    sgn = -sgn;
    N[i] = sgn * lagrange_weight * PolynomialDeCasteljau(weights, n, t);
  }
  return N;
}


double barycentric(vector<double> values, vector<double> weights, int n, double t){
  double N = 0;
  double D = 0;
  for (int i=0; i<=n; i++) {
    double r = t-cos(Pi*i/n);
    if (abs(r) < 1.0e-14)
      return values[i];
    r = weights[i]/r;
    N += r * values[i];
    D += r;
  }
  return N/D;  
}