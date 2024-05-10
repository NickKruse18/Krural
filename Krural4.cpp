#include <Rcpp.h>
using namespace Rcpp;




int KN4Z[100] = {};
double KN4Z2[100] = {};
double KN4x[1000] = {};



// [[Rcpp::export]]
void CloseCenters(const NumericVector B, const int k, const int n, const int c){
  double d = 0;  double M = 0;  int I = 0;  double dd = 0;  const int N = n * k;  int i;  int j;
  for(i = 0; i < c; i++){
    d = std::abs(B[k*i]-KN4x[0]);
    for(j = 1; j < k; j++){
      dd = std::abs(B[k*i+j]-KN4x[j]);
      if(d<dd){ d = dd; }
    }
    KN4Z[i] = i;  KN4Z2[i] = d;
    if(d>M){ M = d;  I = i; }
  }
  i = c*k;
  while(i < N){
    d = std::abs(B[i]-KN4x[0]);
    for(j = 1; j < k; j++){
      i++;  dd = std::abs(B[i]-KN4x[j]);
      if(d<dd){ d = dd; }
    }
    if(d<M){
      KN4Z[I] = i/k;  KN4Z2[I] = d;  M = d;
      for(j = 0; j < c; j++){ if(KN4Z2[j]>M){ M = KN4Z2[j];  I = j; } }
    }
    i++;
  }
}

// [[Rcpp::export]]
NumericVector Krural4Fitted(const NumericVector B, const NumericVector b, const NumericVector X, const int k, const int n, const int m, const int c){
  NumericVector z (n);  double nn;
  int i;  int j;  int I = 0;  double d;
  for(i = 0; i < n; i++){
    for(j = 0; j < k; j++){ KN4x[j] = X[I];  I++; }
    CloseCenters(B,k,m,c);
    nn = 0;  j = 0;  d = 0;
    while(j < c){
      if(KN4Z2[j]==0){ break; }
      d += b[KN4Z[j]]/KN4Z2[j];  nn += 1/KN4Z2[j];  j++; }
    if(j < c){
      nn = 1;  d = b[KN4Z[j]];  j++;
      while(j < c){ if(KN4Z2[j]==0){ d += b[KN4Z[j]];  nn++; };  j++; }
    }
    z[i] = d/nn;
  }
  return z;
}





