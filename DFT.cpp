#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <cstdlib>

using namespace std;


vector<complex<double>> dft(const vector<complex<double>> & x) {
    const int N = x.size();

    vector<complex<double>> X_N(N);

    for(size_t k = 0; k <= N - 1; k++)
        for(size_t n = 0; n <= N - 1; n++)
            X_N[k] = X_N[k] + exp(complex<double>( 0, -2 * M_PI * k * n/ N))* x[n];

    return X_N;

}

vector<complex<double>> randomVector (int n, int upperBound) {
  vector<complex<double>> vec (n);
  for (size_t i = 0; i < vec.size(); i++) {
    vec[i] = complex<double>(random () % upperBound,0);
  }
  return vec;
}

int main(int argc, char *argv[]) {
    // vector<complex<double>> x {complex<double>(2,0),complex<double>(31,0),
    // complex<double>(1,0),complex<double>(0,0),complex<double>(3,0),complex<double>(20,0),
    // complex<double>(4,0),complex<double>(2,0),complex<double>(2,0),complex<double>(31,0),
    // complex<double>(1,0),complex<double>(0,0),complex<double>(3,0),complex<double>(20,0),
    // complex<double>(4,0),complex<double>(2,0)};

    const size_t signal_length = (argc > 1 && atoi(argv[1]) > 0) ? pow(4,atoi(argv[1])) : 4;

    vector<complex<double>> x(signal_length);
    
    srand(40);

    x = randomVector(signal_length, 100);

    vector<complex<double>> X_N(x.size());
    
    const size_t number_of_experiments = 1000;
    ofstream outdata;
    outdata.open("dft_"+to_string((argc > 1 && atoi(argv[1]) > 0) ? (size_t) pow(4,atoi(argv[1])) : 4)+"_"+to_string(number_of_experiments)+".dat");
    double total_duration = 0;
    for(size_t experiment = 0; experiment < number_of_experiments; experiment++ ){
        std::clock_t start;
        double duration;
        start = std::clock();
        X_N = dft(x);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        total_duration += duration;
    }

    outdata << total_duration/number_of_experiments;
    outdata.close();

    // for(auto el : X_N)
    //     std::cout << el << std::endl;

    return 0;
}
