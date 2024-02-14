#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <ctime>

using namespace NTL;

ZZ_pX randomPoly(long n) {

    ZZ_pX a = random_ZZ_pX(n);
    return a;
}

int main() {
    clock_t start, end;
    start = clock();
    long n = 1 << 15;
    ZZ q = NextPrime(power_ZZ(2, 1000));
    ZZ_p::init(q);
    ZZ_pX a, b;
    end = clock();
    std::cout << "Time for setup of polynomial Ring Z mod (~2^1000) of degree 2^15: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;    
    start = clock();
    a = random_ZZ_pX(n);
    end = clock();
    std::cout << "Time for generation of a random element: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
    // std::cout << a << std::endl;
    return 0;
}
