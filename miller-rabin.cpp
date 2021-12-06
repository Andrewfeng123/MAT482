#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>
#include <chrono>

/* Goal: Given a positive integer, test whether it is a prime
*/

typedef boost::multiprecision::uint128_t Integer;
typedef boost::random::uniform_int_distribution<Integer> uniform_int_distrib;
using namespace std::chrono;

Integer gcd(Integer a, Integer b) {
    Integer c;
    if (b < a) {
        c = b;
        b = a;
        a = c;
    }
    while (a != 0) {
        c = a;
        a = b % a;
        b = c;
    }
    return b;
}

// a to the power of b by repeated squaring, arithmetic is done mod c
// if c is unspecified, arithmetic overflow mod 2^128 as limited by Integer
Integer poww(Integer a, Integer b, Integer c = 0) {
    if (b == 0)
        return 1;
    else if (b == 1)
        return c == 0 ? a : a % c;
    else if (b % 2 == 1)
        return c == 0 ? (a * poww((a*a), (b-1)/2)) : ((a % c) * poww(((a*a)%c), (b-1)/2, c)) % c;
    else
        return c == 0 ? poww((a*a), b/2) : poww(((a*a)%c), b/2, c) % c;
}

int k = 20;                 // Here, we fix the number of iterations of MR to 20
                            // So the error rate should be at most 1/(2^40)

bool miller_rabin(Integer n) {
    int d = 0, c = 0; n = n - 1;
    while (n % 2 == 0) {
        n = n/2; d++;
    }
    n = n * poww(2, d) + 1;

    boost::random::mt11213b gen;
    uniform_int_distrib udis(1,n-1);
    while (c < k) {
        Integer b = udis(gen);
        Integer bb = poww(b, n-1, n);
        if (bb != 1 || (gcd(bb - 1, n) > 1 && gcd(bb-1, n) < n))
            return false;
        bb = poww(b, (n-1) / poww(2,d), n);
        for (int j = 0; j < d; j ++) {
            if (gcd(bb - 1, n) > 1 && gcd(bb - 1, n) < n)
                return false;
            bb = (bb * bb) % n;         // hope that bb*bb < 2^128.... (should be good if n < 2^64)
        }
        c++;
    }
    return true;
}

bool naive(Integer n) {
    for (Integer i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            return false;
        }
    }
    return true;
}

int main() {
    Integer n = (poww(2,61)-1) * (poww(2,31)-1);
    std::cout << "Testing whether " << n << " is prime..." << std::endl;
    auto start = high_resolution_clock::now();
    naive(n);
    auto end = high_resolution_clock::now();
    std::cout << "Time taken for naive algorithm: " << duration_cast<milliseconds>(end - start).count() << " milliseconds" << std::endl;

    start = high_resolution_clock::now();
    miller_rabin(n);
    end = high_resolution_clock::now();
    std::cout << "Time taken for miller-rabin algorithm: " << duration_cast<milliseconds>(end - start).count() << " milliseconds" << std::endl;
    getchar();
}