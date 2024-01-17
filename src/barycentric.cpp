#include "barycentric.h"

vector<mpreal> get_homogeneous_values(vector<mpreal> values, vector<mpreal> weights, int n){
  vector<mpreal> Values;
  mpreal t;
  Values.resize(n+1);

  for (int i=0; i<=n; i++){
    Values[i] = values[i] * weights[i];
  }

  vector<mpreal> N;
  N.resize(n+1);
  for (int i=0; i<=n; i++){
    t = cos(Pi*i/n);
    N[i] = PolynomialDeCasteljau(Values, n, t);
  }
  return N;
}

vector<mpreal> get_barycentric_weights(vector<mpreal> weights, int n){
  vector<mpreal> N;
  mpreal t;
  mpreal lagrange_weight;
  int sgn = 1;
  
  N.resize(n+1);
  for (int i=0; i<=n; i++){
    t = cos(Pi*i/n);
    // t = i/n;
    if (i==0 || i==n){
      lagrange_weight = 0.5;
    }else{
      lagrange_weight = 1;
    }
    N[i] = sgn * lagrange_weight * PolynomialDeCasteljau(weights, n, t);
    sgn = -sgn;
  }
  return N;
}

vector<mpreal> get_barycentric_data(vector<mpreal> values, vector<mpreal> weights, int n){
  int sgn = 1;
  mpreal lagrange_weight;
  for (int i=0; i<=n; i++){
    if (i==0 || i==n){
      lagrange_weight = 0.5;
    }else{
      lagrange_weight = 1;
    }
    values[i] = sgn * lagrange_weight * values[i] / weights[i];
    sgn = -sgn;
  }
  return values;
}

mpreal barycentric(vector<mpreal> values, vector<mpreal> weights, int n, mpreal t){
  mpreal N = 0;
  mpreal D = 0;
  for (int i=0; i<=n; i++) {
    mpreal r = t-cos(Pi*i/n);
    if (abs(r) < 1.0e-14)
      return values[i];
    r = weights[i]/r;
    N += r * values[i];
    D += r;
  }
  return N/D;  
}