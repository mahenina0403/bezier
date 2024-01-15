#include "linearGeometric.h"

double linearGeometric(vector<double> values, vector<double> weights, int n, double t)
{
	double h = 1.0;
	double u = 1.0 - t;
	unsigned int n1 = n + 1;
	double R = values[n];
	if(t <= 0.5) {
		u = t / u;
		for(int k=1; k<=n; k++) {
			h = h * u * (n1-k) * weights[n-k];
			h = h / (k * weights[n-k+1] + h);
			double h1 = 1 - h;
			R = h1 * R + h * values[n-k];
		}
	}
	else {
		u = u / t;
		for(int k = 1; k <= n; k++) {
			h = h * (n1-k) * weights[n-k];
			h = h / (k * u * weights[n-k+1] + h);
			double h1 = 1 - h;
			R = h1 * R + h * values[n-k];
		}
	}
	return R;
}