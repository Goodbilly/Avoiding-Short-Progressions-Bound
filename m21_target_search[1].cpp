#include <bits/stdc++.h>
using namespace std;
const int M=21;
double approx(int p,double u,int nb){double worst=0;double x[21];for(int t=0;t<nb;t++){double beta=p*(t+.371)/nb;for(int k=0;k<M;k++){x[k]=fmod(-u*k*k-beta*k,p);if(x[k]<0)x[k]+=p;}sort(x,x+M);double g=p-x[20]+x[0];for(int k=1;k<M;k++)g=max(g,x[k]-x[k-1]);worst=max(worst,g);}return worst;}
int main(){vector<pair<int,int>> qs={{1,4},{1,9},{2,15},{3,4}};for(auto [a,b]:qs){double best=2;tuple<int,int,int,int,double> rec;for(int p=5;p<=85;p++){int s=(p-1)/4;for(int h=1;h<=b*16*5;h++){double R=(double)h/16;double rr=fmod(2*R,p);if(rr<2*s-1e-8||rr>p-2*s+1e-8)continue;double U=(double)a*R/b;double v=approx(p,U,48)/s;if(v<best){best=v;rec={p,s,h,16,v};}if(v<1.000001){cout<<a<<'/'<<b<<" candidate p="<<p<<" s="<<s<<" R="<<h<<"/16 approximate="<<v<<'\n';}}}auto [p,s,h,d,v]=rec;cerr<<a<<'/'<<b<<" best p="<<p<<" s="<<s<<" R="<<h<<'/'<<d<<" sampled="<<v<<'\n';}}
