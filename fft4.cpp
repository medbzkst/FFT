#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <cstdlib>

using namespace std;


vector<complex<double>> fft4(const vector<complex<double>> & x) {
    const int N = x.size();
    if (N <= 1)
        return vector<complex<double>> {x[0]};
    else {
        vector<complex<double>> x_B0(N/4);
        vector<complex<double>> x_B1(N/4);
        vector<complex<double>> x_B2(N/4);
        vector<complex<double>> x_B3(N/4);

        for(size_t k = 0; k <= N/4-1; k++) {
            x_B0[k] = x[4*k];
            x_B1[k]  = x[4*k+1];
            x_B2[k]  = x[4*k+2];
            x_B3[k]  = x[4*k+3];
        }

        vector<complex<double>> X_B0_Ndiv4 = fft4(x_B0);
        vector<complex<double>> X_B1_Ndiv4 = fft4(x_B1);
        vector<complex<double>> X_B2_Ndiv4 = fft4(x_B2);
        vector<complex<double>> X_B3_Ndiv4 = fft4(x_B3);

        vector<complex<double>> X_N(N);

        for(size_t k = 0; k <= N/4 - 1; k++) {

            const complex<double> WkN = exp(complex<double>( 0, -2 * M_PI * k / N));
            const complex<double> W2kN = exp(complex<double>( 0, -2 * M_PI * 2*k / N));
            const complex<double> W3kN = exp(complex<double>( 0, -2 * M_PI * 3*k / N));

            X_N[k]        = X_B0_Ndiv4[k]    +                      WkN*X_B1_Ndiv4[k]  + W2kN*X_B2_Ndiv4[k]  +                      W3kN*X_B3_Ndiv4[k];
            X_N[k+ N/4]   = X_B0_Ndiv4[k]    - complex<double>(0,1)*WkN*X_B1_Ndiv4[k]  - W2kN*X_B2_Ndiv4[k]  + complex<double>(0,1)*W3kN*X_B3_Ndiv4[k];
            X_N[k+ N/2]   = X_B0_Ndiv4[k]    -                      WkN*X_B1_Ndiv4[k]  + W2kN*X_B2_Ndiv4[k]  -                      W3kN*X_B3_Ndiv4[k];
            X_N[k+ 3*N/4] = X_B0_Ndiv4[k]    + complex<double>(0,1)*WkN*X_B1_Ndiv4[k]  - W2kN*X_B2_Ndiv4[k]  - complex<double>(0,1)*W3kN*X_B3_Ndiv4[k];
        }

        return X_N;

    }

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
    outdata.open("fft4_"+to_string((argc > 1 && atoi(argv[1]) > 0) ? (size_t) pow(4,atoi(argv[1])) : 4)+"_"+to_string(number_of_experiments)+".dat");
    double total_duration = 0;
    for(size_t experiment = 0; experiment < number_of_experiments; experiment++ ){
        std::clock_t start;
        double duration;
        start = std::clock();
        X_N = fft4(x);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        total_duration += duration;
    }

    outdata << total_duration/number_of_experiments;
    outdata.close();
    // for(auto el : X_N)
    //     std::cout << el << std::endl;

    return 0;
}