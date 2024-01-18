#include "main.h"

void compare_runtime(vector<mpreal> f, vector<mpreal> beta, int n, int sample){
	mpreal t;
	int i;
	clock_t startTime;
	clock_t endTime;

	ofstream sampletime("../data/result.csv", ios::app);

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalDeCasteljau(f,beta,n,t);
	}
	endTime = clock();
	// cout << "Method: RationalDeCasteljau" << endl;
	cout << "RDC: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		FarinRationalDeCasteljau(f,beta,n,t);
	}
	endTime = clock();
	// cout << "Method: FarinRationalDeCasteljau" << endl;
	cout << "FDC: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalVS(f,beta,n,t);
	}
	endTime = clock();
	// cout << "Method: RationalVS" << endl;
	cout << "RVS: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalHornBez(f,beta,n,t);
	}
	endTime = clock();
	
	// cout << "Method: RationalHornBez" << endl;
	cout << "RHB: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		RationalLader(f,beta,n,t);
	}
	endTime = clock();
	// cout << "Method: RationalLader" << endl;
	cout << "RLD: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		linearGeometric(f,beta,n,t);
	}
	endTime = clock();
	// cout << "Method: linearGeometric" << endl;
	cout << "LTG: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	auto bar_values = get_homogeneous_values(f,beta,n);
	auto bar_weights = get_barycentric_weights(beta,n);
	bar_values = get_barycentric_data(bar_values, bar_weights,n);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		barycentric(bar_values,bar_weights,n,t);
	}
	endTime = clock();
	// cout << "Method: barycentric" << endl;
	cout << "BAR: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << ",";

	startTime = clock();
	auto W = convert_to_wang_ball(f,beta,n);
	for(i=0; i<=sample; i++){
		t = 1.0*i/sample;
		rationalWangBall(W,n,t);
	}
	endTime = clock();
	// cout << "Method: rationalWangBall" << endl;
	cout << "RWB: Time elapsed is " << (endTime-startTime) / (double) CLOCKS_PER_SEC << " s" << endl;
	sampletime << (endTime-startTime) / (double) CLOCKS_PER_SEC << endl;

	sampletime.close();
}

void compare_precision(vector<mpreal> f, vector<mpreal> beta, int n, mpreal t){
	mpreal res;
	res = RationalDeCasteljau(f,beta,n,t);
	cout << "FDC: " << res - FarinRationalDeCasteljau(f,beta,n,t) << endl;
	cout << "RVS: " << res - RationalVS(f,beta,n,t) << endl;
	cout << "RHB: " << res - RationalHornBez(f,beta,n,t) << endl;
	cout << "RLD: " << res - RationalLader(f,beta,n,t) << endl;
	cout << "LTG: " << res - linearGeometric(f,beta,n,t) << endl;

	auto bar_values = get_homogeneous_values(f,beta,n);
	auto bar_weights = get_barycentric_weights(beta,n);
	bar_values = get_barycentric_data(bar_values, bar_weights,n);
	cout << "BAR: " << res - barycentric(bar_values,bar_weights,n,t) << endl;

	auto W = convert_to_wang_ball(f,beta,n);
	cout << "RWB: " << res - rationalWangBall(W,n,t) << endl;

}

int main(int argc, char *argv[]) {

	mpreal::set_default_prec(1024);
	
	int n = 3;
	int sample = 1000;
	
	mpreal t = 0.5;

	vector<mpreal> f;
	vector<mpreal> beta;
	f.resize(n+1);
	beta.resize(n+1);

	f[0] = 0.1;
	f[1] = 0.00000002;
	f[2] = 0.0000001;
	f[3] = 0.1;

	beta[0] = 1;
	beta[1] = 0.1;
	beta[2] = 0.3;
	beta[3] = 1;

	if (argc == 2 && strcmp(argv[1],"-t") == 0){
		int Samples[6] = {1000, 10000, 20000, 30000, 40000, 50000};
		ofstream sampletime("../data/result.csv");
		sampletime.close();
		for (auto s : Samples){
			compare_runtime(f,beta,n,s);	
		}
	}else if (argc == 3 && strcmp(argv[1],"-t") == 0){
		sample = atoi(argv[2]);
		compare_runtime(f,beta,n,sample);
	}else if (argc == 2 && strcmp(argv[1],"-p") == 0){
		compare_precision(f,beta,n,t);
	}else if (argc == 3 && strcmp(argv[1],"-p") == 0){
		t = atof(argv[2]);
		compare_precision(f,beta,n,t);
	}else{
		cerr << "Command error" << endl;
	}
	return 0;
}