#include "main.h"

void compare_runtime(const vector<vec2> f, const vector<double> beta, int n, int sample, ofstream& sampletime){
	double t;
	int i;
	clock_t startTime;
	clock_t endTime;

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);
	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalDeCasteljau(f,beta,n,t);
	}
	endTime = clock();
	cout << "RDC: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		FarinRationalDeCasteljau(f,beta,n,t);
	}
	endTime = clock();
	cout << "FDC: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";


	startTime = clock();
	for(i=0;i<=n;i++){
		g[i] = f[i];
		alpha[i] = beta[i];
	}
	gen_VS_data(&g, &alpha, n);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalVS2(g,alpha,n,t);
	}
	endTime = clock();
	cout << "RVS: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	// startTime = clock();
	// for(i=0; i<=sample; i++){
	// 	t = 1.0*i/sample;
	// 	RationalVS(f,beta,n,t);
	// }
	// endTime = clock();
	// cout << "VSo: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	// sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalHornBez(f,beta,n,t);
	}
	endTime = clock();

	cout << "RHB: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		linearGeometric(f,beta,n,t);
	}
	endTime = clock();
	cout << "LTG: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	vector<double> T;
	// startTime = clock();
	// T = compute_nodes(n,CHEBYSHEV);
	// get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	// for(i=0; i<=sample; i++){
	// 	t = 1.0*i/sample;
	// 	barycentric(g,alpha,T,n,t);
	// }
	// endTime = clock();
	// // cout << "Method: barycentric" << endl;
	// cout << "CH1: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	// sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	
	for(i=0; i<=sample/2; i++){
		t = 1.0*i/sample;
		barycentric_2(g,alpha,T,n,t);
	}
	endTime = clock();
	cout << "BCH: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	// T = compute_nodes(n,UNIFORM);
	// get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	// for(i=0; i<=sample; i++){
	// 	t = 1.0*i/sample;
	// 	barycentric(g,alpha,T,n,t);
	// }
	// endTime = clock();
	// cout << "UN1: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	// sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	for(i=0; i<=sample/2; i++){
		t = 1.0*i/sample;
		barycentric_2(g,alpha,T,n,t);
	}
	endTime = clock();
	// cout << "Method: barycentric" << endl;
	cout << "BUN: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	startTime = clock();
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		rationalWangBall(WangBallPoints,WangBallWeights,t);
	}
	endTime = clock();
	cout << "RWB: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	vector<complex<double>> x(n+1);
	vector<complex<double>> y(n+1);
	vector<complex<double>> z(n+1);
	vector<complex<double>> r(n+1);

	startTime = clock();
	for (int i=0; i<= n; i++){
		x[i] = beta[i]*f[i].x();
		y[i] = beta[i]*f[i].y();
		z[i] = beta[i];
	}

	x = toFourier(x,n);
	y = toFourier(y,n);
	z = toFourier(z,n);

	r = roots_of_unity(n);

	for(i=0; i<=sample/2; i++){
		t = 1.0*i/sample;
		BernsteinFourrier_2(x,y,z,r,n,t);
	}
	endTime = clock();
	cout << "BFF: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC<< " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << endl;
}


void compare_values(const vector<vec2> f,const vector<double> beta, int n, double t){

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	for(int i=0;i<=n;i++){
		g[i] = f[i];
		alpha[i] = beta[i];
	}
	gen_VS_data(&g, &alpha, n);

	// cout << "Homogeneous-Binomial data" << endl;
	// for(int i=0; i<=n; i++){
	// 	cout << g[i] << " " << alpha[i] << endl;
	// }
	
	cout << "RDC: " << RationalDeCasteljau(f,beta,n,t) << endl;
	cout << "FDC: " << FarinRationalDeCasteljau(f,beta,n,t) << endl;
	cout << "RVS: " << RationalVS2(g,alpha,n,t) << endl;
	cout << "RHB: " << RationalHornBez(g,alpha,n,t) << endl;
	cout << "LGM: " << linearGeometric(f,beta,n,t) << endl;
		

	vector<double> T;
	
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	auto u = barycentric_2(g,alpha,T,n,t);
	cout << "CHE: " << u[0] << endl;
	
	// cout << "barycentric CHEBYSHEV data" << endl;
	// for(int i=0; i<=n; i++){
	// 	cout << g[i] << " " << alpha[i] << endl;
	// }

	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	u = barycentric_2(g,alpha,T,n,t);
	cout << "UNI: " << u[0] << endl;

	// cout << "barycentric UNIFORM data" << endl;
	// for(int i=0; i<=n; i++){
	// 	cout << g[i] << " " << alpha[i] << endl;
	// }

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	cout << "RWB: " << rationalWangBall(WangBallPoints,WangBallWeights,t) << endl;

	vector<complex<double>> x(n+1);
	vector<complex<double>> y(n+1);
	vector<complex<double>> z(n+1);
	vector<complex<double>> r(n+1);

	for (int i=0; i<= n; i++){
		x[i] = beta[i]*f[i].x();
		y[i] = beta[i]*f[i].y();
		z[i] = beta[i];
	}

	x = toFourier(x,n);
	y = toFourier(y,n);
	z = toFourier(z,n);

	r = roots_of_unity(n);

	cout << "BFF: " << BernsteinFourrier(x,y,z,r,n,t) << endl;
}


int main(int argc, char *argv[]) {

	int n;
	int sample = 100;
	string filename = "result.csv";

	if (argc == 2)
		sample = atoi(argv[1]);

	double t = 0.8;

	vector<vec2> f;
	vector<double> beta;
	
	n = 20;
	f.resize(n+1);
	beta.resize(n+1);
	for(int i=0; i<=n; i++){
		f[i] = vec2(i*100);
		beta[i] = (i%2)+1;
	}

	vector<int> S{100,1000,10000,100000,500000,1000000};

	if (argc==1){
		ofstream sampletime("");
		compare_runtime(f,beta,n,sample,sampletime);
	}
	else if (argc==2 && strcmp(argv[1],"-f")==0){
		ofstream sampletime(filename);
		sampletime.close();
		for (const auto& s: S){
			ofstream sampletime(filename, ios::app);
			compare_runtime(f,beta,n,s,sampletime);
			sampletime.close();
		}
	}
	else if (argc==3 && strcmp(argv[1],"-v")==0){
		t = atof(argv[2]);
		compare_values(f,beta,n,t);
	}
	else if (argc==2){
		sample = atoi(argv[1]);
		ofstream sampletime("");
		compare_runtime(f,beta,n,sample,sampletime);
	}
	
	return 0;
}