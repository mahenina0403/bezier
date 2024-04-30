// SPDX-License-Identifier: MIT
// Copyright (c) [2024] [Fuda Chiara, Andriamahenina Ramanantoanina]


#include "main.h"

void to_of(vec2 p, ofstream& of, double eol=0){
	of << p.x() << "'" << p.y();
	if (eol) of << endl;
	else of << ", ";
}

void to_of(mpreal x, mpreal y, ofstream& of, double eol=0){
	of << x << "'" << y;
	if (eol) of << endl;
	else of << ",";
}

void compare_runtime(const vector<vec2> f, const vector<double> beta, int n, int sample, ofstream& sampletime, double show=false){
	double t;
	int i;
	clock_t startTime;
	clock_t endTime;
	int population = 1000;
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
	if (show) cout << "RDC: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "FDC: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "RVS: Time elapsed is " << runtime << " s" << endl;
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

	if (show)  cout << "RHB: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "LTG: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "BCH: Time elapsed is " << runtime << " s" << endl;
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

	if (show)  cout << "BUN: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "RWB: Time elapsed is " << runtime << " s" << endl;
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
	if (show)  cout << "BFF: Time elapsed is " << runtime << " s" << endl;
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
	// cout << "RWB: " << rationalWangBall(WangBallPoints,WangBallWeights,t) << endl;

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


void compare_error(const vector<vec2> f,const vector<double> beta, int n, double t, ofstream& sampletime, double show=false){

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	gen_VS_data(f, beta, &g, &alpha, n);

	vector<double> fx(n+1);
	vector<double> fy(n+1);
	for(int i=0; i<=n; i++){
		fx[i] = f[i].x();
		fy[i] = f[i].y();
	}

	mpreal exact_valuex = RationalDeCasteljau(fx,beta,n,t);
	mpreal exact_valuey = RationalDeCasteljau(fy,beta,n,t);
	
	mpreal exact_x = abs(exact_valuex);
	mpreal exact_y = abs(exact_valuey);
	
	vec2 rdc = RationalDeCasteljau(f,beta,n,t);
	vec2 fdc = FarinRationalDeCasteljau(f,beta,n,t);
	vec2 rvs = RationalVS2(g,alpha,n,t);
	vec2 rhb = RationalHornBez(g,alpha,n,t);
	vec2 ltg = linearGeometric(f,beta,n,t);

	vector<double> T;
	
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	vec2 che = barycentric(g,alpha,T,n,t);
	
	
	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	vec2 uni = barycentric(g,alpha,T,n,t);
	
	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	vec2 rwb = rationalWangBall(WangBallPoints,WangBallWeights,t);

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
	vec2 rbf = BernsteinFourrier(x,y,z,r,n,t);
	
	mpreal X;
	mpreal Y;

	// cout << exact_y << endl;
    X = abs(exact_valuex - (mpreal)rdc.x())/exact_x; Y= abs(exact_valuey - (mpreal)rdc.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)fdc.x())/exact_x; Y= abs(exact_valuey - (mpreal)fdc.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)rvs.x())/exact_x; Y= abs(exact_valuey - (mpreal)rvs.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)rhb.x())/exact_x; Y= abs(exact_valuey - (mpreal)rhb.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)ltg.x())/exact_x; Y= abs(exact_valuey - (mpreal)ltg.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)che.x())/exact_x; Y= abs(exact_valuey - (mpreal)che.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)uni.x())/exact_x; Y= abs(exact_valuey - (mpreal)uni.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)rwb.x())/exact_x; Y= abs(exact_valuey - (mpreal)rwb.y())/exact_y;
    to_of(X,Y,sampletime);
    X = abs(exact_valuex - (mpreal)rbf.x())/exact_x; Y= abs(exact_valuey - (mpreal)rbf.y())/exact_y;
	to_of(X,Y,sampletime,true);
    
    if (show){
    	cout.precision(5);
    	X = abs(exact_valuex - (mpreal)rdc.x())/exact_x; Y= abs(exact_valuey - (mpreal)rdc.y())/exact_y;
    	cout << setw(6) << left << "rdc:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)fdc.x())/exact_x; Y= abs(exact_valuey - (mpreal)fdc.y())/exact_y;
	    cout << setw(6) << left << "fdc:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)rvs.x())/exact_x; Y= abs(exact_valuey - (mpreal)rvs.y())/exact_y;
	    cout << setw(6) << left << "rvs:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)rhb.x())/exact_x; Y= abs(exact_valuey - (mpreal)rhb.y())/exact_y;
	    cout << setw(6) << left << "rhb:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)ltg.x())/exact_x; Y= abs(exact_valuey - (mpreal)ltg.y())/exact_y;
	    cout << setw(6) << left << "ltg:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)che.x())/exact_x; Y= abs(exact_valuey - (mpreal)che.y())/exact_y;
	    cout << setw(6) << left << "che:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)uni.x())/exact_x; Y= abs(exact_valuey - (mpreal)uni.y())/exact_y;
	    cout << setw(6) << left << "uni:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)rwb.x())/exact_x; Y= abs(exact_valuey - (mpreal)rwb.y())/exact_y;
	    cout << setw(6) << left << "rwb:";
	    cout << setw(13) << X << setw(20) << Y << endl;
	    X = abs(exact_valuex - (mpreal)rbf.x())/exact_x; Y= abs(exact_valuey - (mpreal)rbf.y())/exact_y;
		cout << setw(6) << left << "rbf:";
	    cout << setw(13) << X << setw(20) << Y << endl;
    }   
}

void rcond(const vector<vec2> f,const vector<double> beta, int n, double t, ofstream& sampletime){
	vector<vec2> g(n+1);
	vector<double> alpha(n+1);

	vec2 u = rcond_rdc(f,beta,n,t);
	to_of(u,sampletime);

	vector<double> T;
	
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	u = rcond_barycentric(g,alpha,T,n,t);
	to_of(u,sampletime);

	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	u = rcond_barycentric(g,alpha,T,n,t);
	to_of(u,sampletime);

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	u = rcond_WangBall(WangBallPoints,WangBallWeights,t);
	to_of(u,sampletime,true);

}


int main(int argc, char *argv[]) {
	int my_mpreal_precision = 1024;
    mpreal::set_default_prec(my_mpreal_precision);

	int n;
	int sample = 100;
	string filename = "result.csv";

	if (argc == 2)
		sample = atoi(argv[1]);

	double pi = Pi();

	double t = 0.5;

	vector<vec2> f;
	vector<double> beta;


	n = 50;
	f.resize(n+1);
	beta.resize(n+1);
	for(int i=0; i<=n; i++){
		// example of unstable for UNI, RWB, RBF
		f[i] = vec2(1, sin(i*pi/(n))+1);
		if (i<10 || i>n-10) f[i].x() = pow(10,6);
		beta[i] = (i%2)+1;
	}

	n = 4;
	f.resize(n+1);
	beta.resize(n+1);
	for(int i=0; i<=n; i++){
		// example of unstable for UNI, RWB, RBF
		beta[i] = 1;
	}
	f[0] = vec2(10, -100);
	f[1] = vec2(20, 200);
	f[2] = vec2(30, -200);
	f[3] = vec2(40, 101);
	f[4] = vec2(50, 101);

	// vector<int> S{1, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750};
	vector<int> S{1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750};
	vector<int> degrees{3,5,7,10,13,15,17,20};

	if (argc==1){
		n = int(rand()*10)+3;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			f[i] = vec2(rand());
			beta[i] = rand()+1;
		}
		ofstream sampletime("");
		compare_runtime(f,beta,n,1000,sampletime,true);
	}
	else if (argc==2 && strcmp(argv[1],"-rd")==0){
		filename = "runtime-degree.csv";
		ofstream sampletime(filename);
		sampletime.close();
		for (int s=3; s<=20; s++){
			n = s;
			f.resize(n+1);
			beta.resize(n+1);
			for(int i=0; i<=n; i++){
				// f[i] = vec2(i*100+1, sin(i*pi/(n))+1);
				f[i] = 100*i*vec2(1)+vec2(1);
				beta[i] = (i%2)+1;
			}
			ofstream sampletime(filename, ios::app);
			compare_runtime(f,beta,n,2500,sampletime);
			sampletime.close();
		}
	}
	else if (argc==2 && strcmp(argv[1],"-rs")==0){
		n = 20;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			// example of unstable for UNI, RWB, RBF
			f[i] = 100*i*vec2(1)+vec2(1);
			beta[i] = (i%2)+1;
		}
		filename = "runtime-sample.csv";
		ofstream rsd20(filename);
		rsd20.close();
		for (const auto& s: S){
			ofstream rsd20(filename, ios::app);
			compare_runtime(f,beta,n,s,rsd20);
			rsd20.close();
		}
	}
	else if (argc==2 && strcmp(argv[1],"-e")==0){
		ofstream x("relative-error.csv");
		x.close();
	
		for(int i=0; i<= 1000; i++){
			t = i/1000.0;
			ofstream x("relative-error.csv", ios::app); 
			compare_error(f,beta,n,t,x);
			x.close();
		}

	}
	else if (argc==3 && strcmp(argv[1],"-e")==0){
		t = atof(argv[2]);
		ofstream x("");
		x.close();
	
		compare_error(f,beta,n,t,x,true);
		x.close();
	}
	else if (argc==2 && strcmp(argv[1],"-c")==0){
		ofstream x("rcond.csv");
		x.close();
	
		for(int i=0; i<= 1000; i++){
			t = i/1000.0;
			ofstream x("rcond.csv", ios::app); 
			rcond(f,beta,n,t,x);
			x.close();
		}

	}
	else if (argc==3 && strcmp(argv[1],"-r")==0){
		n = int(rand()*10)+3;
		f.resize(n+1);
		beta.resize(n+1);
		for(int i=0; i<=n; i++){
			f[i] = vec2(rand());
			beta[i] = rand()+1;
		}

		sample = atoi(argv[2]);
		ofstream sampletime("");
		compare_runtime(f,beta,n,sample,sampletime,true);
	}

	return 0;
}