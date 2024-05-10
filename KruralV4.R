




Krural4 = function(X,Y,par=2,C=5){
  S = sum(Y^2)
  k = ncol(X)
  n = nrow(X)
  m = 0
  B = b = c()
  for(i in 1:k){
    
  }
  for(i in 1:nrow(X)){
    x = X[i,];  y = Y[i]
    O = optim(par=c(x,y),fn=function(z){
      r = Krural4Fitted(c(B,z[1:k]),c(b,z[k+1]),X,k,n,m+1,min(C,m+1))
      return(sum((Y-r)^2))
    }, method="BFGS")
    O$value =  O$value + par * (m + 1)
    print(i)
    print(c(O$value,S))
    print(O$par)
    if(O$value < S){ B = c(B,O$par[1:k]);  b = c(b,O$par[k+1]);  S = O$value;  m = m + 1 }
    if(i==10){ break }
  }
  print(Y)
  print(Krural4Fitted(c(B,10,10),c(b,100),X,k,n,m+1,min(C,m+1)))
  print(c(sum((Y-Krural4Fitted(c(B,10,10),c(b,100),X,k,n,m+1,min(C,m+1)))^2)+ par * (m + 1),S))
  return(list(B,b,k))
}

Krural4 = function(X,Y,par=2,start=10,C=5){
  S = sum(Y^2)
  k = ncol(X)
  n = nrow(X)
  m = start
  XX = as.vector(t(X))
  if(start==n){ return(list(XX,Y,k)) }
  set.seed(1)
  I = sample(1:n,start)
  x = c(as.vector(t(X[I,,drop=F])),Y[I])
  print(sum((Y-Krural4Fitted(x[1:(k*start)],x[k*start+1:start],XX,k,n,m,min(C,m)))^2))
  O = optim(par=x,fn=function(z){
    r = Krural4Fitted(z[1:(k*start)],z[k*start+1:start],XX,k,n,m,min(C,m))
    return(sum((Y-r)^2))
  }, method="BFGS")
  print(O$value/n)
  return(list(O$par[1:(k*start)],O$par[k*start+1:start],k))
}

B = Krural4(Data[[1]][1:900,],Data[[2]][1:900],2,900,5)
Data[[2]][1:900]

mean((Data[[2]][1:900]-Krural4Fitted(B[[1]],B[[2]],as.vector(t(Data[[1]][1:900,])),2,900,900,5))^2)
mean((Data[[2]][901:1000]-Krural4Fitted(B[[1]],B[[2]],as.vector(t(Data[[1]][901:1000,])),2,100,900,5))^2)



Data[[2]][1:1000]-Krural4Fitted(as.vector(t(Data[[1]][1:900,])),Data[[2]][1:900],as.vector(t(Data[[1]][1:1000,])),2,1000,800,10)

microbenchmark::microbenchmark(Krural4Predict(list(Data[[1]][1:900,],Data[[2]][1:900]),Data[[1]][1:900,]))

mean((Data[[2]][1:900]-Krural4Predict(list(Data[[1]][1:900,],Data[[2]][1:900]),Data[[1]][1:900,]))^2)

mean((Data[[2]][901:1000]-Krural4Predict(list(Data[[1]][1:900,],Data[[2]][1:900]),Data[[1]][901:1000,]))^2)
mean((Data[[2]][9001:10000]-Krural4Predict(list(Data[[1]][1:9000,],Data[[2]][1:9000]),Data[[1]][9001:10000,]))^2)

Data[[2]][901:1000]-Krural4Predict(list(Data[[1]][1:900,],Data[[2]][1:900]),Data[[1]][901:1000,])
Data[[2]][901:1000]-predict(xgb,Data[[1]][901:1000,])


Data[[1]][999,]

hist(Data[[2]][901:1000]-Krural4Predict(list(Data[[1]][1:900,],Data[[2]][1:900]),Data[[1]][901:1000,]),100)

xgb = xgboost::xgboost(Data[[1]][1:900,],Data[[2]][1:900],nrounds=100)

mean((Data[[2]][901:1000]-predict(xgb,Data[[1]][901:1000,]))^2)


hist(Data[[2]][901:1000]-predict(xgb,Data[[1]][901:1000,]),100)






B = as.vector(t(Data[[1]][1:900,]))
b = as.vector(t(Data[[1]][801:1000,]))

Rcpp::sourceCpp("Krural4.cpp")
mean((Data[[2]][1:1000]-Krural4Fitted(as.vector(t(Data[[1]][1:900,])),Data[[2]][1:900],as.vector(t(Data[[1]][1:1000,])),2,1000,900,10))^2)
microbenchmark::microbenchmark(Krural4Fitted(B,Data[[2]][1:900],b,2,200,900,10),times=10000)
