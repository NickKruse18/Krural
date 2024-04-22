
Rcpp::sourceCpp("Segment.cpp")

TopOptim = function(r,B){
  O = optim(par = B,fn = TopCalcScore, r = r, method = "BFGS")
  return(O)
}

PredictKruralNickwork = function(B,X){
  z = X%*%B[[1]]
  InitializeKrural(X[1:length(X)],nrow(X),ncol(X))
  for(i in 1:nrow(B[[2]])){
    z = z + TopCalc(B[[2]][i,])
  }
  return(z)
}

KruralNickwork = function(X,y,overfitting = 1){
  k = ncol(X);  n = nrow(X)
  InitializeKrural(X[1:length(X)],n,k)
  H = solve(t(X)%*%X,t(X))
  b = H%*%y;  r = y-X%*%b;  s = sum(r^2)
  print(var(r))
  opt = 2*var(y) * overfitting
  I = which.max(abs(r))
  S = numeric(k)
  for(i in 1:k){ S[i] = max(X[,i]) - min(X[,i]) }
  ind = 1;  i = 1
  Tops = matrix(0,1,2+2*k)
  Tops[i,] = c(r[I],r[I],X[I,],r[I]*S/4)
  Z = matrix(0,n,1)
  Z[,i] = TopCalc(Tops[i,]);
  while(T){
    Tops[i,] = TopOptim(r,Tops[i,])$par
    Z[,i] = TopCalc(Tops[i,])
    o = y - Z%*%rep(1,ncol(Z));  b = H%*%o;  o = o-X%*%b
    if(sum(o^2)<s-opt){
      s = sum(o^2);  r = y - X%*%b;  ind = 0
      for(j in 1:nrow(Tops)){
        Tops[j,] = TopOptim(r-Z[,-j,drop=F]%*%rep(1,ncol(Z)-1),Tops[j,])$par
        Z[,j] = TopCalc(Tops[j,])
      }
      i = i + 1
      o = y - Z%*%rep(1,ncol(Z))
      b = H%*%o;  r = o - X%*%b
      Tops = rbind(Tops,0);  Z = cbind(Z,0)
      print(var(r))
    }
    else if(ind == n){ break }
    ind = ind + 1
    I = match(sort(r)[ind],r)
    Tops[i,] = c(r[I],r[I],X[I,],r[I]*S/4)
    Z[,i] = TopCalc(Tops[i,]);
  }
  return(list(b,Tops))
}





































