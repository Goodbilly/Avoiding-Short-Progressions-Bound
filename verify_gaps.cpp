// Recalculate every claimed G_{p,21}(j/10) in the CSVs using exact arithmetic.
#define main original_M22_verifier_main
#include "gap_core.cpp"
#undef main
#include <fstream>
#include <sstream>
#include <string>
int main() {
    int verified=0;
    for (const char *filename : {"odd_moduli_exact_gaps.csv","even_moduli_exact_gaps.csv"}) {
        std::ifstream in(filename);
        if (!in) { std::cerr << "Missing " << filename << "\n"; return 2; }
        std::string line;
        while (std::getline(in,line)) {
            if(line.empty()) continue;
            std::replace(line.begin(),line.end(),',',' ');
            std::istringstream ss(line);
            int m,p,s,j; long long gn,gd;
            if (!(ss>>m>>p>>s>>j>>gn>>gd)) return 3;
            if(m!=21 || gn>s*gd) continue;
            Rational actual=G(p,m,j,10);
            if ((i128)actual.numerator*gd!=(i128)gn*actual.denominator) {
                std::cerr<<"Mismatch: "<<p<<" "<<j<<"/10\n";
                return 1;
            }
            ++verified;
        }
    }
    std::cout<<"Exact gap entries verified: "<<verified<<"\n";
}
