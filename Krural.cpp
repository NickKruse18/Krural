#include <Rcpp.h>
using namespace Rcpp;

double KNX[1000000] = {};
int KNn = 0;
int KNk = 0;

// [[Rcpp::export]]
void InitializeKrural(NumericVector X, int n, int k){
  KNn = n;  KNk = k;  for(int i = 0; i < n*k; i++){ KNX[i] = X[i]; }
}

// [[Rcpp::export]]
NumericVector TopCalc(NumericVector B){
  double a = 0;  int I = 0;
  NumericVector z (KNn);
  for(int i = 0; i < KNn; i++){
    a = B[1];  I = i;
    for(int j = 0; j < KNk; j++){
      a = std::min(a, B[0] - B[2+KNk+j] * std::abs(B[2+j] - KNX[I]));
      if(a <= 0){ a = 0;  break; }
      I += KNn;
    }
    z[i] = a;
  }
  return z;
}

// [[Rcpp::export]]
double TopCalcScore(NumericVector B, NumericVector r){
  double s = 0;  double a = 0;  NumericVector z (KNn);  z = TopCalc(B);
  for(int i = 0; i < KNn; i++){ a = r[i]-z[i];  s += a*a; }
  return s;
}


