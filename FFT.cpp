#include <iostream>
#include <complex>
#include <chrono>
#include <random>
#include <boost/multiprecision/cpp_complex.hpp>
#include <cmath>

/* Goal: Given a list of n=2^m complex numbers representing an n-dimensional 
    complex vector, return the vector in the basis in discrete fourier transform

*/

typedef boost::multiprecision::cpp_complex_quad Complex;

// fast
void FFT(std::vector<Complex>& in, std::vector<Complex>& out);

// elementary
void eFFT(std::vector<Complex>& in, std::vector<Complex>& out);

Complex i = Complex(0.0,1.0);
Complex pi = Complex(boost::math::constants::pi<long double>(), 0);

int main() {
    int N = pow(2, 5); int j = 1, k = 1;    
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 10.0);

    std::vector<std::vector<Complex>> ins(k);
    // initialize k fft inputs
    for (int kk = 0; kk < k; kk++) {
        ins[kk] = std::vector<Complex>(N);
        for (int l = 0; l < N; l++) {
            ins[kk][l] = {dist(mt),dist(mt)};
        }
    }

    int fft_ms = 0, efft_ms = 0;    // store total time taken for FFT and eFFT
    std::vector<Complex> out(N);
    std::vector<Complex> out1(N);

    for (int jj = 0; jj < j; jj++) {
        using namespace std::chrono;

        auto start = high_resolution_clock::now();
        for (int kk = 0; kk < k; kk++) {
            FFT(ins[kk], out);
        }
        auto end = high_resolution_clock::now();
        fft_ms += duration_cast<milliseconds>(end - start).count();

        start = high_resolution_clock::now();
        for (int kk = 0; kk < k; kk++) {
            eFFT(ins[kk], out1);
        }
        end = high_resolution_clock::now();
        efft_ms += duration_cast<milliseconds>(end - start).count();

        // std::cout << std::setprecision(2);
        // for(Complex & c : out)
        //     std::cout << c;
        // std::cout << std::endl;
    }

    std::cout << "Performing random DFT on length " << N << " vectors" << std::endl;
    std::cout << "FFT: " << fft_ms << " milliseconds; " << "elementary DFT: " 
                << efft_ms << " milliseconds" << std::endl;
    
    double diff = 0;
    for (int kk = 0; kk < N; kk++) {
        diff = std::max(diff, norm(out[kk] - out1[kk]).real().convert_to<double>());
    }
    std::cout << "max entry-wise difference: " << diff << std::endl;

    getchar();
    return 0;
}

// FFT
// in: length N array of complex numbers where N is a power of 2
// out: length N array that is the discrete fourier transform of N, preallocated
//      to size N

void FFT(std::vector<Complex>& in, std::vector<Complex>& out) {
    int N = in.size();
    if (N == 1) {
        out[0] = in[0];
    } else {
        std::vector<Complex> P(N/2), Q(N/2), out_P(N/2), out_Q(N/2);
        for (int j = 0; j < N/2; j++) {
            P[j] = in[2*j];
            Q[j] = in[2*j + 1];
        }

        // evaluate P and Q at the (N/2)th roots of unity
        FFT(P, out_P);
        FFT(Q, out_Q);

        // Compute final result
        for (int j = 0; j < N; j++) {
            out[j] = 0.5 * (out_P[j % (N/2)] 
                        + exp(-2*pi*i*j / N) * out_Q[j % (N/2)]);
        }
    }
}

// Elementary School Algorithm
// input and output same as above
void eFFT(std::vector<Complex>& in, std::vector<Complex>& out) {
    int N = in.size();
    for(int j = 0; j < N; j++) {
        out[j] = 0;
        for(int k = 0; k < N; k++) {
            out[j] += in[k] * exp(-2*pi*i*j*k / N);
        }
        out[j] = out[j] / N;
    }
}
