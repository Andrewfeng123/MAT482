# MAT482
Selected algorithms from MAT482 at the University of Toronto, math topics course titled Algorithms In Number Theory and Algebra. The code is written in C++, using Clang as the compiler, Boost::multiprecision library for arbitrary precision arithmetic (set to a limited precision in the code, but could be better if desired). The command for compilation is written in FFT.txt and only needs slight modification to work for miller-rabin.

## Fast Fourier Tranform
One of the coolest and useful algorithm. Can compute the discrete Fourier Transform of a length n list of complex number in O(n log n) time!

## Miller-Rabin Primality Test
Another very useful algorithm. This algorithm is randomized and correct with high probability (as high as you want!). If the input is a prime, then the algorithm definitely says "prime"; if the input is composite, then the algorithm can prove it in polynomial time with high probability! For more detailed information, check out Miller-Rabin.pdf.
