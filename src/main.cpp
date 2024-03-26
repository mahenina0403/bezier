#include "main.h"

void compare_runtime(const vector<vec2> f, const vector<double> beta, int n, int sample, ofstream& sampletime){
	double t;
	int i;
	clock_t startTime;
	clock_t endTime;

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

	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalVS(f,beta,n,t);
	}
	endTime = clock();
	cout << "RVS: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

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

	vector<vec2> g(n+1);
	vector<double> alpha(n+1);
	
	vector<double> T;
	startTime = clock();
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		barycentric(g,alpha,T,n,t);
	}
	endTime = clock();
	// cout << "Method: barycentric" << endl;
	cout << "BAR CHE: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		barycentric(g,alpha,T,n,t);
	}
	endTime = clock();
	cout << "BAR UNI: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
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
	// cout << "Method: rationalWangBall" << endl;
	cout << "RWB: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << endl;
}


void compare_values(const vector<vec2> f,const vector<double> beta, int n, double t){

	cout << "RDC: " << RationalDeCasteljau(f,beta,n,t) << endl;
	cout << "FDC: " << FarinRationalDeCasteljau(f,beta,n,t) << endl;
	cout << "RVS: " << RationalVS(f,beta,n,t) << endl;
	cout << "RHB: " << RationalHornBez(f,beta,n,t) << endl;
	cout << "LGM: " << linearGeometric(f,beta,n,t) << endl;
		
	vector<vec2> g(n+1);
	vector<double> alpha(n+1);
	vector<double> T;
	
	T = compute_nodes(n,CHEBYSHEV);
	get_data(f,beta,T,n,&g,&alpha,CHEBYSHEV);
	cout << "CHE: " << barycentric(g,alpha,T,n,t) << endl;
	
	T = compute_nodes(n,UNIFORM);
	get_data(f,beta,T,n,&g,&alpha,UNIFORM);
	cout << "UNI: " << barycentric(g,alpha,T,n,t) << endl;

	vector<vec2> WangBallPoints(n+1);
	vector<double> WangBallWeights(n+1);
	
	convert_to_wang_ball(f,beta,&WangBallPoints,&WangBallWeights,n);
	cout << "RWB: " << rationalWangBall(WangBallPoints,WangBallWeights,t) << endl;
}


int main(int argc, char *argv[]) {

	int n;
	int sample = 100;
	string filename = "../data/sample_result.csv";

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
		beta[i] = 1;
	}

	vector<int> S{100,1000,10000,100000,500000,1000000};


	if (argc==2){
		ofstream sampletime("");
		compare_runtime(f,beta,n,sample,sampletime);
		return 0;
	}

	ofstream sampletime(filename);
	sampletime.close();
	for (const auto& s: S){
		ofstream sampletime(filename, ios::app);
		compare_runtime(f,beta,n,s,sampletime);
		sampletime.close();
	}

	// if (argc==2)
	// 	t = atof(argv[1]);
	// compare_values(f,beta,n,t);
	
	return 0;
}