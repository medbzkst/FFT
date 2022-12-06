# Fast Fourier Transform
This repository contains Fast Fourier Transform algorithm implementation in C++ using DFT (is the baseline, but not a FAST algorithm), Cooley-Tukey, and Radix-4 algorithm.

The files `DFT.cpp`, `fft2.cpp`, `fft4.cpp` are the source code of DFT, Cooley-Tukey, and Radix-4 respectively. They are all self contained for simplicity. Each of which has the function that implements the corresponding algorithm as the file name is. They all contain a function `randomVector` that generates random double values in a vector of `std::complex` placed in the real part. They all use the same seed for fair comparison. The function `main` in each file prepares the vector that represents the signal to transform to the frequency space, applies the corresponding algorithm for a hardwired number of times (hardwired to 1000) and outputs the mean execution time in a file.

The experiment flow is about to run the executable with an argument `K` where `K` will be the exponent of 4 to define the number of elements that constitute the signal, i.e. if `./fft2 5` is executed, the Fourier Transform for a signal of length `4^K` which is `4^5` is computed using the Cooley-Tukey algorithm. Notably, the values of that constitutes the complex vector of the signal do not exceed 100. This value is hardwired.

## Compilation

`gcc` compiler is used with the following arguments `-lstdc++ -lm -ldl`, i.e. to compile Cooley-Tukey in file `fft2.cpp`, the command to execute is `gcc fft2.cpp -lstdc++ -lm -ldl -o fft2`.

## Output files

Each output file will have the name `ALGO_sizeOfSignal_NumberOfExperiments.dat` where `ALGO` is the name of the executed function, `sizeOfSignal` is the number of elements consituting the signal vector, and `NumberOfExperiments` is the hardwired 1000 times of performing the experiment, i.e. `fft4_64_1000.dt` contains the required _mean_ execution time to run Radix-4 algorithm 1000 times on a signal of 64 elements.
