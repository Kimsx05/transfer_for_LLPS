#include <Rcpp.h>
using namespace Rcpp;
// Exact type-7 quartiles used by CellChat::triMean. No approximation or sampling.
// [[Rcpp::export]]
NumericMatrix grouped_trimean_cpp(NumericMatrix x, IntegerVector group, int ng) {
 const int nr=x.nrow(),nc=x.ncol();
 std::vector<std::vector<int>> indices(ng);
 for(int i=0;i<nr;++i){if(group[i]<1||group[i]>ng)stop("Invalid group");indices[group[i]-1].push_back(i);}
 NumericMatrix out(ng,nc);
 for(int j=0;j<nc;++j){
  if(j%100==0)checkUserInterrupt();
  for(int g=0;g<ng;++g){
   int n=indices[g].size();if(n==0)stop("Empty factor level");std::vector<double> v(n);
   for(int k=0;k<n;++k)v[k]=x(indices[g][k],j);
   std::sort(v.begin(),v.end());
   double q[3];double pp[3]={0.25,0.5,0.75};
   for(int k=0;k<3;++k){double pos=(n-1)*pp[k];int lo=(int)std::floor(pos);int hi=(int)std::ceil(pos);double h=pos-lo;q[k]=(lo==hi)?v[lo]:((1-h)*v[lo]+h*v[hi]);}
   out(g,j)=(q[0]+q[1]+q[1]+q[2])/4.0;
  }
 }
 return out;
}
