#include "main.h"

void compare_runtime(const vector<vec2> f, const vector<double> beta, int n, int sample, ofstream& sampletime){
	double t;
	int i;
	clock_t startTime;
	clock_t endTime;
	int population = 100;
	double runtime;
	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			RationalDeCasteljau(f,beta,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "RDC: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			FarinRationalDeCasteljau(f,beta,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "FDC: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		gen_VS_data(f, beta, &g, &alpha, n);
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			RationalVS2(g,alpha,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "RVS: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		gen_VS_data(f, beta, &g, &alpha, n);
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			RationalHornBez(f,beta,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;

	cout << "RHB: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			linearGeometric(f,beta,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "LTG: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	vector<double> T;

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		T = compute_nodes(n,CHEBYSHEV);
		get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
		
		for(i=0; i<=sample/2; i++){
			t = 1.0*i/sample;
			barycentric_2(g,alpha,T,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "BCH: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		T = compute_nodes(n,UNIFORM);
		get_data(f,beta,T,n,&g,&alpha,UNIFORM);
		for(i=0; i<=sample/2; i++){
			t = 1.0*i/sample;
			barycentric_2(g,alpha,T,n,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;

	cout << "BUN: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	runtime = 0;
	for(int l=0; l<= population; l++){
		startTime = clock();
		convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
		for(i=0; i<=sample; i++){
			t = 1.0*i/sample;
			rationalWangBall(WangBallPoints,WangBallWeights,t);
		}
		endTime = clock();
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "RWB: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << ",";

	vector<complex<double>> x(n+1);
	vector<complex<double>> y(n+1);
	vector<complex<double>> z(n+1);
	vector<complex<double>> r(n+1);

	runtime = 0;
	for(int l=0; l<= population; l++){
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
		runtime += (endTime-startTime) / (double) CLOCKS_PER_SEC;
	}
	runtime /= population;
	cout << "BFF: Time elapsed is " << runtime << " s" << endl;
	sampletime << runtime << endl;
}


void compare_values(const vector<vec2> f,const vector<double> beta, int n, double t){

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	gen_VS_data(f, beta, &g, &alpha, n);

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

	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	u = barycentric_2(g,alpha,T,n,t);
	cout << "UNI: " << u[0] << endl;

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
	u = BernsteinFourrier_2(x,y,z,r,n,t);
	cout << "BFF: " << u[0] << endl;
}


void compare_error(const vector<vec2> f,const vector<double> beta, int n, double t, ofstream& sampletime){

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	gen_VS_data(f, beta, &g, &alpha, n);

	vector<double> fx(n+1);
	for(int i=0; i<=n; i++)
		fx[i] = f[i].x();

	mpreal exact_value = RationalDeCasteljau(fx,beta,n,t);
	// cout << RationalDeCasteljau(fx,beta,n,t) - FarinRationalDeCasteljau(fx,beta,n,t) << endl;

	double rdc = RationalDeCasteljau(f,beta,n,t).x();
	double fdc = FarinRationalDeCasteljau(f,beta,n,t).x();
	double rvs = RationalVS2(g,alpha,n,t).x();
	double rhb = RationalHornBez(g,alpha,n,t).x();
	double ltg = linearGeometric(f,beta,n,t).x();

	vector<double> T;
	
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	auto u = barycentric_2(g,alpha,T,n,t);
	double che = u[0].x();
	
	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	u = barycentric_2(g,alpha,T,n,t);
	double uni = u[0].x();

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	double rwb = rationalWangBall(WangBallPoints,WangBallWeights,t).x();

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
	u = BernsteinFourrier_2(x,y,z,r,n,t);
	double bff = u[0].x();
	

	cout.precision(5);
    cout << setw(50) << left << "DeCasteljau:";
    cout << (exact_value - (mpreal)rdc)/exact_value << endl;
    cout << setw(50) << "FarinDeCasteljau: ";
    cout << (exact_value - (mpreal)fdc)/exact_value << endl;
    cout << setw(50) << "VS:";
    cout << (exact_value - (mpreal)rvs)/exact_value << endl;
    cout << setw(50) << "HornBez:";
    cout << (exact_value - (mpreal)rhb)/exact_value << endl;
    cout << setw(50) << "linearGeometric: ";
    cout << (exact_value - (mpreal)ltg)/exact_value << endl;
    cout << setw(50) << "barycentric (UNIFORM):";
    cout << (exact_value - (mpreal)uni)/exact_value << endl;
    cout << setw(50) << "barycentric (CHEBYSHEV):";
    cout << (exact_value - (mpreal)che)/exact_value << endl;
    cout << setw(50) << "Wang-Ball:";
    cout << (exact_value - (mpreal)rwb)/exact_value << endl;
    cout << setw(50) << "Bernstein-Fourier:";
    cout << (exact_value - (mpreal)bff)/exact_value << endl;


    sampletime << (exact_value - (mpreal)rdc)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)fdc)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)rvs)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)rhb)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)ltg)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)che)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)uni)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)rwb)/exact_value << ", ";
    sampletime << (exact_value - (mpreal)bff)/exact_value << endl;
}

int main(int argc, char *argv[]) {
	int my_mpreal_precision = 64;
    mpreal::set_default_prec(my_mpreal_precision);

	int n;
	int sample = 100;
	string filename = "result.csv";

	if (argc == 2)
		sample = atoi(argv[1]);

	double t = 0.3;

	vector<vec2> f;
	vector<double> beta;
	
	n = 3;
	f.resize(n+1);
	beta.resize(n+1);
	for(int i=0; i<=n; i++){
		// f[i] = vec2(i*100)+vec2(1,0);
		f[i] = vec2(rand()%100+1);
		beta[i] = (i%2)+1;
	}

	vector<int> S{250, 500, 1000,2000,5000,10000};
	vector<int> degrees{3,5,7,10,15,20};

	if (argc==1){
		ofstream sampletime("");
		compare_runtime(f,beta,n,1000,sampletime);
	}
	else if (argc==2 && strcmp(argv[1],"-rd")==0){
		filename = "runtime-degree.csv";
		ofstream sampletime(filename);
		sampletime.close();
		for (const auto& s: degrees){
			n = s;
			f.resize(n+1);
			beta.resize(n+1);
			for(int i=0; i<=n; i++){
				f[i] = vec2(i*100)+vec2(1,0);
				beta[i] = (i%2)+1;
			}
			ofstream sampletime(filename, ios::app);
			compare_runtime(f,beta,n,5000,sampletime);
			sampletime.close();
		}
	}
	else if (argc==2 && strcmp(argv[1],"-rs")==0){
		n = 3;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			// f[i] = vec2(i*100)+vec2(1,0);
			f[i] = vec2(rand()%100+1);
			beta[i] = (i%2)+1;
		}

		filename = "runtime-sample-deg-3.csv";
		ofstream rsd3(filename);
		rsd3.close();
		for (const auto& s: S){
			ofstream rsd3(filename, ios::app);
			compare_runtime(f,beta,n,s,rsd3);
			rsd3.close();
		}

		n = 20;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			// f[i] = vec2(i*100)+vec2(1,0);
			f[i] = vec2(rand()%100+1);
			beta[i] = (i%2)+1;
		}
	
		filename = "runtime-sample-deg-20.csv";
		ofstream rsd20(filename);
		rsd20.close();
		for (const auto& s: S){
			ofstream rsd20(filename, ios::app);
			compare_runtime(f,beta,n,s,rsd20);
			rsd20.close();
		}
	}
	else if (argc==3 && strcmp(argv[1],"-e")==0){
		t = atof(argv[2]);
		ofstream sampletime("");
		compare_error(f,beta,n,t,sampletime);
	}
	else if (argc==2 && strcmp(argv[1],"-e")==0){
		n = 3;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			f[i] = vec2(i*100)+vec2(1,0);
			// f[i] = vec2(rand()%100+1);
			beta[i] = (i%2)+1;
		}

		filename = "relative-error-deg-3.csv";
		ofstream ed3(filename);
		ed3.close();
		for(int i=0; i<= 250; i++){
			t = i/250.0;
			ofstream ed3(filename, ios::app);
			compare_error(f,beta,n,t,ed3);
			ed3.close();
		}

		n = 20;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			f[i] = vec2(i*100)+vec2(1,0);
			// f[i] = vec2(rand()%100+1);
			beta[i] = (i%2)+1;
		}

		filename = "relative-error-deg-20.csv";
		ofstream ed20(filename);
		ed20.close();
		for(int i=0; i<= 250; i++){
			t = i/250.0;
			ofstream ed20(filename, ios::app);
			compare_error(f,beta,n,t,ed20);
			ed20.close();
		}

	}
	else if (argc==3 && strcmp(argv[1],"-r")==0){
		sample = atoi(argv[2]);
		ofstream sampletime("");
		compare_runtime(f,beta,n,sample,sampletime);
	}

	return 0;
}