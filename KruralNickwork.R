

Rcpp::sourceCpp("Krural.cpp")

InitializeKrural(Xiris[1:600],150,4)

microbenchmark::microbenchmark(TopCalc(c(9,9,Xiris[1,],c(3,3,3,3))))


TopCalcScore(Iris[,5],c(9,9,Xiris[1,],c(3,3,3,3)))

TopCalc(c(9,9,Xiris[1,],c(3,3,3,3)))
TopCalc(c(-9,-9,Xiris[1,],c(3,3,3,3)))

TopCalc(c(-100,-100,-10,10,5,5))
TopCalc(c(100,100,-10,10,10,10))

Xiris

microbenchmark::microbenchmark(sum((Iris[,5] - TopCalc(c(9,9,Xiris[1,],c(3,3,3,3))))^2))
microbenchmark::microbenchmark(TopCalcScore(Iris[,5],c(9,9,Xiris[1,],c(3,3,3,3))))

microbenchmark::microbenchmark(TopOptim(Iris[,6],c(1,1,Xiris[51,],c(1,1,1,1))))


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

TopOptim(Iris[,6],Xiris,c(1,0,Xiris[51,],c(1,1,1,1)))


Iris[,6]-round(TopCalc(TopOptim(Iris[,6],Xiris,c(1,0,Xiris[51,],c(1,1,1,1)))$par,Xiris))

plot(Xiris[,c(1,2)],col=1+Iris[,6])

KruralNickwork = function(X,y,overfitting = 2){
  k = ncol(X);  n = nrow(X)
  InitializeKrural(X[1:length(X)],n,k)
  H = solve(t(X)%*%X,t(X))
  b = H%*%y;  r = y-X%*%b;  s = sum(r^2)
  rr = var(r)
  print(var(r))
  opt = overfitting
  I = which.max(abs(r))
  S = numeric(k)
  for(i in 1:k){ S[i] = max(X[,i]) - min(X[,i]);  if(S[i]==0){ S[i] = 1 } }
  ind = 1;  i = 1
  Tops = matrix(0,1,2+2*k)
  Tops[i,] = c(r[I],r[I],X[I,],abs(r[I])/S*4)
  Z = matrix(0,n,1)
  Z[,i] = TopCalc(Tops[i,]);
  while(T){
    Tops[i,] = TopOptim(r,Tops[i,])$par
    Z[,i] = TopCalc(Tops[i,])
    o = y - Z%*%rep(1,ncol(Z));  b = H%*%o;  o = o-X%*%b
    if((mean(r^2)-mean(o^2))/rr>opt){
      s = sum(o^2);  r = y - X%*%b;  ind = 0
      j = 1;
      while(j <= nrow(Tops)){
        o = r - Z[,-j,drop=F]%*%rep(1,ncol(Z)-1)
        oo = mean(o^2)
        Tops[j,] = TopOptim(o,Tops[j,])$par
        Z[,j] = TopCalc(Tops[j,])
        o = o - Z[,j]
        if((oo-mean(o^2))/rr<opt){ Z = Z[,-j];  Tops = Tops[-j,];  i = i - 1;  j = j - 1;  print(j) }
        j = j + 1;
      }
      i = i + 1
      o = y - Z%*%rep(1,ncol(Z))
      b = H%*%o;  r = o - X%*%b
      Tops = rbind(Tops,0);  Z = cbind(Z,0)
      print(c(mean(r^2),nrow(Tops)))
    }
    if(ind == n){ break }
    R = r[order(abs(r),decreasing = T)]
    while(T){
      ind = ind + 1
      I = match(R[ind],r)
      t = c(r[I],r[I],X[I,],abs(r[I])/S*4)
      z = TopCalc(t)
      break
    }
    Tops[i,] = t
    Z[,i] = z
  }
  return(list(b,Tops))
}



B = KruralNickwork(Xiris,Iris[,6],10)

round(PredictKruralNickwork(B,Xiris))-Iris[,5]

B
Rcpp::sourceCpp("Krural.cpp")




































