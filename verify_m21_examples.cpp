#include <bits/stdc++.h>
using namespace std;
#define main original_four_main
#include "verify_four.cpp"
#undef main
int main(){
 struct E{int a,b,p,s; long long rn,rd;};
 vector<E> examples={{1,4,30,7,7,1},{1,9,37,9,19,2},{2,15,19,4,23,1},{3,4,45,11,67,2},{1,4,64,15,15,1},{1,9,53,13,213,16},{2,15,71,17,147,8},{3,4,24,5,19,1}};
 for(auto e:examples){auto g=G(e.p,21,e.a*e.rn,e.b*e.rd);cout<<e.a<<'/'<<e.b<<" p="<<e.p<<" s="<<e.s<<" R="<<e.rn<<'/'<<e.rd<<" G="<<g.numerator<<'/'<<g.denominator<<" pass="<<((__int128)g.numerator<= (__int128)e.s*g.denominator)<<'\n';}
}
