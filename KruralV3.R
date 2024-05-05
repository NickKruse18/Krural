


Bifurcate = function(X,y,SM=NULL,A=FALSE){
  if(A){ K = kmeans(X,2) }
  else{ K = kmeans(y,2) }
  if(is.null(SM)){ SM = matrix(0,2,ncol(X));  for(i in 1:ncol(X)){ SM[,i] = c(min(X[,i]),max(X[,i])) } }
  S = matrix(0,4,ncol(X))
  N = numeric(4)
  xI = matrix(T,nrow(X),2)
  for(I in 1:2){
    XK = X[K$cluster==I,,drop=F]
    for(i in 1:ncol(XK)){ S[2*I-2+1:2,i] = c(min(XK[,i]),max(XK[,i])) }
    
    for(i in 1:ncol(XK)){ xI[,I] = xI[,I]&(X[,i]>=S[2*I-1,i])&(X[,i]<=S[2*I,i]) }
    N[2*I-2+1:2] = c(sum((y[xI[,I]]-mean(y[xI[,I]]))^2),sum((y[!xI[,I]]-mean(y[!xI[,I]]))^2))
  }
  I = which.min(c(N[1]+N[2],N[3]+N[4]))
  S = S[2*I-2+1:2,]
  xI = xI[,I]
  xi = sum(xI)
  for(i in 1:ncol(X)){
    xII = rep(T,nrow(X))
    for(j in 1:ncol(X)){ if(j==i){ next };  xII = xII&(X[,j]>=S[1,j])&(X[,j]<=S[2,j]) }
    s = SM[1,i] - S[1,i]
    xIII = xII&(X[,i]<=S[2,i])
    while(sum(xIII&(X[,i]>=S[1,i]+s))!=xi){ s = s/2 }
    S[1,i] = S[1,i] + s
    s = SM[2,i] - S[2,i]
    xIII = xII&(X[,i]>=S[1,i])
    while(sum(xIII&(X[,i]<=S[2,i]+s))!=xi){ s = s/2 }
    S[2,i] = S[2,i] + s
  }
  #plot(X[,1:2],col=1+xI,xlim=c(0,10),ylim=c(0,10))
  #lines(c(SM[1,1],SM[2,1],SM[2,1],SM[1,1],SM[1,1]),c(SM[1,2],SM[1,2],SM[2,2],SM[2,2],SM[1,2]),col=1)
  #lines(c(S[1,1],S[2,1],S[2,1],S[1,1],S[1,1]),c(S[1,2],S[1,2],S[2,2],S[2,2],S[1,2]),col=2)
  S = rbind(S,SM)
  return(list(S,N[2*I-2+1:2],xI))
}

Bifurcate(Xiris,Iris[,5])
Bifurcate(Data[[1]][1:900,],Data[[2]][1:900])
plot(Xiris[,c(3,4)],col=1+Iris[,5])


Krural3 = function(X,y,par=2,S=NULL,n=NULL){
  if(is.null(n)){ n = sum((y-mean(y))^2) }
  if(is.null(S)){ S = matrix(0,2,ncol(X));  for(i in 1:ncol(X)){ S[,i] = c(min(X[,i]),max(X[,i])) } }
  if(n<par){ return(list(S,mean(y),n)) }
  if(length(y)<3){ return(list(S,mean(y),n)) }
  B = Bifurcate(X,y,S)
  if(sum(B[[1]][1:2,]==S)==0){ B = Bifurcate(X,y,S,T) }
  N = B[[2]]
  if(sum(N)+par>n){ return(list(S,mean(y),n)) }
  S = B[[1]]
  xI = B[[3]]
  B1 = Krural3(X[xI,],y[xI],par,S[1:2,],N[1])
  B2 = Krural3(X[!xI,],y[!xI],par,S[3:4,],N[2])
  return(list(rbind(B1[[1]],B2[[1]]),c(B1[[2]],B2[[2]]),B1[[3]]+B2[[3]]))
}

PredictKrural3 = function(B,X){
  P = numeric(nrow(X))
  S = B[[1]]
  for(i in 1:nrow(X)){
    x = X[i,]
    n = nrow(S)
    j = 1
    while(j < n){
      W = T
      for(k in 1:ncol(X)){
        if(x[k]<S[j,k]){ W = F;  break }
        if(x[k]>S[j+1,k]){ W = F;  break }
      }
      if(W){ P[i] = B[[2]][(j+1)/2];  break }
      j = j + 2
    }
  }
  return(P)
}

VisualizeKrural = function(B,n){
  Im = matrix(0,n,n)
  k = nrow(B[[1]])/2
  x = c(B[[1]][2*k-1,1],B[[1]][2*k,1])
  y = c(B[[1]][2*k-1,2],B[[1]][2*k,2])
  x = (x[2]-x[1])*(0:(n-1))/(n-1) + x[1]
  y = (y[2]-y[1])*(0:(n-1))/(n-1) + y[1]
  print(x)
  print(y)
  for(i1 in 1:n){
    for(i2 in 1:n){
      Im[i1,i2] = PredictKrural3(B,matrix(c(x[i1],y[i2]),1,2))
    }
  }
  image(Im,col = hcl.colors(100, "YlOrRd", rev = TRUE))
}

VisualizeXGB = function(xgb,B,n){
  Im = matrix(0,n,n)
  k = nrow(B[[1]])/2
  x = c(B[[1]][2*k-1,1],B[[1]][2*k,1])
  y = c(B[[1]][2*k-1,2],B[[1]][2*k,2])
  x = (x[2]-x[1])*(0:(n-1))/(n-1) + x[1]
  y = (y[2]-y[1])*(0:(n-1))/(n-1) + y[1]
  print(x)
  print(y)
  for(i1 in 1:n){
    for(i2 in 1:n){
      Im[i1,i2] = predict(xgb,matrix(c(x[i1],y[i2]),1,2))
    }
  }
  image(Im,col = hcl.colors(100, "YlOrRd", rev = TRUE))
}
?image
B = Krural3(Xiris,Iris[,5],1)
Iris[,5]-PredictKrural3(B,Xiris)


VisualizeKrural(B,101)
VisualizeXGB(xgb,B,101)

Data[[1]]

B = Krural3(Data[[1]][1:9000,],Data[[2]][1:9000],0.01)
B[[3]]/9000
length(B[[2]])

mean((Data[[2]][9001:10000]-PredictKrural3(B,Data[[1]][9001:10000,]))^2)
B

hist(Data[[2]][901:1000]-PredictKrural3(B,Data[[1]][901:1000,]))

cbind(Data[[1]][901:1000,],PredictKrural3(B,Data[[1]][901:1000,]),Data[[2]][901:1000])

xgb = xgboost::xgboost(Data[[1]][1:9000,],Data[[2]][1:9000],nrounds=200)

mean((Data[[2]][9001:10000]-predict(xgb,Data[[1]][9001:10000,]))^2)
mean((Data[[2]][1:9000]-predict(xgb,Data[[1]][1:9000,]))^2)

hist(Data[[2]][9001:10000]-predict(xgb,Data[[1]][9001:10000,]),100)
hist(Data[[2]][9001:10000]-PredictKrural3(B,Data[[1]][9001:10000,]),100)
