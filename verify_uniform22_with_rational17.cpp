// Exact verification for the uniform 22-term Euclidean Ramsey theorem.
//
// Compiling and running:
//   g++ -O2 -std=c++17 verify_uniform22_with_rational17.cpp -o verify
//   ./verify
//
// This program verifies the values of G_{p,22} used in the proof,
// checks G_{33,22}(u) <= 8 throughout [9/25,16/25], and checks
// the 187 remaining reduced rational parameters a/b, and also verifies
// the 20-term result and the rational-interval 17-term proposition.
// All comparisons use integer arithmetic; no floating point is used.
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

// The Lipschitz constant for M=22 is 110. Check [L/d,H/d]
// by evaluating G at the midpoint; bisect only if necessary.
constexpr i64 C = 110;
int evaluations=0, intervals=0;
bool verify_interval(int L, int H, int d) {
    ++evaluations;
    Rational g = G(33, 22, (i64)L+H, 2LL*d);
    i64 margin = 8*g.denominator - g.numerator;
    if (margin < 0) return false;
    if ((i128)2*d*margin >= (i128)C*(H-L)*g.denominator) {
        ++intervals;
        return true;
    }
    // On the indicated input, the process terminates after 63 calls.
    return verify_interval(2*L, L+H, 2*d) &&
           verify_interval(L+H, 2*H, 2*d);
}

i64 floor_div(i64 a, i64 b) {
    assert(b>0);
    i64 q=a/b, r=a%b;
    return q-(r<0);
}
i64 ceil_div(i64 a, i64 b) { return -floor_div(-a,b); }

// Exact interval-intersection test for p=33,s=8,
// red interval [8,17/2], blue interval [9/25,16/25], period 33/2.
bool intersects_large_interval(int a, int b) {
    return ceil_div(18LL*b-425LL*a,825LL) <=
           floor_div(32LL*b-400LL*a,825LL);
}

struct TestValue {int p,s,n,d; i64 gap_n,gap_d;};

// Eleven further rational values used to cover the finite cases.
const vector<TestValue> test_values = {
    {17,4,1,1,4,1},
    {17,4,2,1,4,1},
    {17,4,1,2,4,1},
    {17,4,2,3,133,39},
    {17,4,3,5,39,10},
    {29,6,7,1,6,1},
    {29,6,14,1,17,3},
    {29,6,1,4,23,4},
    {29,6,157,11,731,132},
    {33,8,15,4,15,2},
    {33,8,23,4,8,1}
};

// Exact test for the blue interval around u=n/d with radius
// delta=(s-G(u))/110. The red interval is [s,p/2-s] modulo p/2.
bool intersects_test_interval(int a, int b, const TestValue &z) {
    i64 margin = (i64)z.s*z.gap_d-z.gap_n;
    assert(margin >= 0);
    i64 common = (i64)z.d*C*z.gap_d;
    i64 numerator_left = 2LL*b*((i64)z.n*C*z.gap_d-(i64)z.d*margin)
                      -(i64)a*(z.p-2*z.s)*common;
    i64 numerator_right = 2LL*b*((i64)z.n*C*z.gap_d+(i64)z.d*margin)
                       -(i64)a*(2*z.s)*common;
    i64 denominator = (i64)z.p*common;
    return ceil_div(numerator_left,denominator) <=
           floor_div(numerator_right,denominator);
}

// Verify the three rational values that certify the blue interval
// [7/16,179/400] for M=17, p=37, s=9. The Lipschitz constant is 64.
bool rational_leq(Rational x, Rational y) {
    return (i128)x.numerator*y.denominator <=
           (i128)y.numerator*x.denominator;
}
struct Interval { Rational left, right; };
Interval interval17(i64 n, i64 d, i64 gn, i64 gd) {
    Rational g=G(37,17,n,d);
    assert((i128)g.numerator*gd == (i128)gn*g.denominator);
    i64 margin=9*g.denominator-g.numerator;
    assert(margin>0);
    i64 den=d*64*g.denominator;
    return {reduce(n*64*g.denominator-d*margin,den),
            reduce(n*64*g.denominator+d*margin,den)};
}
void verify_seventeen_interval() {
    Interval I1=interval17(11,25,221,25);
    Interval I2=interval17(111,250,222,25);
    Interval I3=interval17(67,150,134,15);
    assert(rational_leq(Rational{7,16}, I1.left));
    assert(rational_leq(I1.left,I2.left));
    assert(rational_leq(I2.left,I1.right));
    assert(rational_leq(I2.left,I3.left));
    assert(rational_leq(I3.left,I2.right));
    assert(rational_leq(Rational{179,400},I3.right));
    cout << "17-term blue interval [7/16,179/400] verified"
         << " (three exact gap values and overlapping neighborhoods)\n";
}

int main() {
    // Verify every rational value used in the proof.
    for (const TestValue &z : test_values) {
        Rational g = G(z.p,22,z.n,z.d);
        assert((i128)g.numerator*z.gap_d ==
               (i128)z.gap_n*g.denominator);
    }
    Rational g = G(29,22,1,2);
    assert(g.numerator == 17 && g.denominator == 3);
    cout << "Exact gap values verified: " << test_values.size()+1 << "\n";

    // Recover the 20-term result and verify the entire 17-term interval.
    Rational g20 = G(29,20,7,1);
    assert((i128)10*g20.numerator == (i128)67*g20.denominator);
    verify_seventeen_interval();
    cout << "The 20-term gap value verified\n";

    assert(verify_interval(9,16,25));
    assert(evaluations == 63 && intervals == 32);
    cout << "Interval [9/25,16/25] verified: " << evaluations
         << " gap evaluations and " << intervals << " subintervals\n";

    int count=0, direct=0, further=0, missing=0;
    cout << "a | candidates | equation (11) | further coefficients | uncovered\n";
    for (int a=1; a<=5; ++a) {
        int row_total=0, row_direct=0, row_further=0;
        for (int b=1; 25LL*a+14LL*b<825; ++b) {
            if (gcd(a,b)!=1) continue;
            ++count; ++row_total;
            if (intersects_large_interval(a,b)) {
                ++direct; ++row_direct;
                continue;
            }
            bool found=false;
            for (const TestValue &z:test_values) {
                if (intersects_test_interval(a,b,z)) { found=true;break; }
            }
            if (found) { ++further; ++row_further; }
            else { ++missing;cerr << "Not verified: " << a << '/' << b << '\n'; }
        }
        cout << a << " | " << row_total << " | " << row_direct
             << " | " << row_further << " | "
             << (row_total-row_direct-row_further) << "\n";
    }
    // The counts above are computed from the inequalities, not hard-coded.
    // Only the mathematically necessary completeness condition is checked.
    if (missing != 0) {
        cerr << "ERROR: some rational ratios were not certified.\n";
        return 1;
    }
    cout << "Reduced fractions: " << count << "\n";
    cout << "Covered by [9/25,16/25]: " << direct << "\n";
    cout << "Covered by eleven further values: " << further << "\n";
    cout << "Unverified fractions: " << missing << "\n";
    cout << "All checks passed.\n";
}
