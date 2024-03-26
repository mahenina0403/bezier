#include "linearGeometric.h"

vec2 linearGeometric(const vector<vec2> values, const vector<double> weights, int n, double t)
{
	double h = 1.0;
	double u = 1.0 - t;
	unsigned int n1 = n + 1;
	vec2 R = values[0];
	if(t <= 0.5) {
		u = t / u;
		for(int k=1; k<=n; k++) {
			h = h * u * (n1-k) * weights[k];
			h = h / (k * weights[k-1] + h);
			double h1 = 1 - h;
			R = h1 * R + h * values[k];
		}
	}
	else {
		u = u / t;
		for(int k = 1; k <= n; k++) {
			h = h * (n1-k) * weights[k];
			h = h / (k * u * weights[k-1] + h);
			double h1 = 1 - h;
			R = h1 * R + h * values[k];
		}
	}
	return R;
}