// Exact independent verification of four new M=21 certificates.
// Compile: g++ -O2 -std=c++17 verify_four.cpp -o verify_four
// Run: ./verify_four
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using i64 = long long;
using i128 = __int128_t;

struct Rational { i64 numerator, denominator; };
Rational reduce(i64 n, i64 d) {
    assert(d > 0);
    i64 h = gcd(n >= 0 ? n : -n, d);
    return {n/h, d/h};
}

// The maximum circular gap for -u k^2 - beta k, k=0,...,M-1.
// It suffices to evaluate beta at pairwise coincidences.
Rational G(int p, int M, i64 numerator, i64 denominator) {
    assert(p >= 2 && M >= 2 && denominator > 0);
    Rational u = reduce(numerator, denominator);
    i64 U = u.numerator, d = u.denominator;
    U %= (i64)p*d;
    if (U < 0) U += (i64)p*d;
    Rational answer{0,1};
    vector<i64> points(M);
    for (int h=1; h<M; ++h) {
        for (int i=0; i+h<M; ++i) {
            for (int m=0; m<h; ++m) {
                i64 modulus = (i64)p*d*h;
                i64 beta_numerator = (i64)p*m*d - U*h*(2*i+h);
                for (int k=0; k<M; ++k) {
                    i64 residue = (i64)(-((i128)U*h*k*k +
                                           (i128)beta_numerator*k) % modulus);
                    if (residue < 0) residue += modulus;
                    points[k] = residue;
                }
                sort(points.begin(), points.end());
                i64 largest = modulus - points.back() + points.front();
                for (int k=1; k<M; ++k)
                    largest = max(largest, points[k]-points[k-1]);
                i64 gap_denominator = d*h;
                if ((i128)largest*answer.denominator >
                    (i128)answer.numerator*gap_denominator)
                    answer = {largest, gap_denominator};
            }
        }
    }
    return reduce(answer.numerator, answer.denominator);
}


struct Certificate {int a,b,p,s; long long rn,rd,expected_n,expected_d;};
int main(){
  const vector<Certificate> certificates={
    {2,23,5,1,1,1,203,207},
    {1,8,9,2,23,2,2,1},
    {1,10,33,8,115,2,8,1},
    {2,11,5,1,27,2,75,77}
  };
  for(auto z:certificates){
    assert(4*z.s<z.p);
    // R modulo T=p/2 belongs to [s,T-s].
    // Equivalently, 2R modulo p lies in [2s,p-2s].
    long long residue=(2*z.rn)%(z.p*z.rd);
    assert(residue>=2LL*z.s*z.rd && residue<=1LL*(z.p-2*z.s)*z.rd);
    Rational g=G(z.p,21,1LL*z.a*z.rn,1LL*z.b*z.rd);
    assert((i128)g.numerator*z.expected_d==(i128)z.expected_n*g.denominator);
    assert((i128)g.numerator<= (i128)z.s*g.denominator);
    cout<<z.a<<"/"<<z.b<<" p="<<z.p<<" s="<<z.s<<" R="<<z.rn<<"/"<<z.rd
        <<" G="<<g.numerator<<"/"<<g.denominator<<" CERTIFIED\n";
  }
}
