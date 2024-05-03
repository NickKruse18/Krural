


Bifurcate = function(X,y,SM=NULL){
  K = kmeans(y,2)
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
  S = rbind(S[2*I-2+1:2,],SM)
  return(list(S,N[2*I-2+1:2],xI[,I]))
}

Bifurcate(Xiris,Iris[,5])

Krural3 = function(X,y,par=2,S=NULL,n=NULL){
  if(is.null(n)){ n = sum((y-mean(y))^2) }
  if(is.null(S)){ S = matrix(0,2,ncol(X));  for(i in 1:ncol(X)){ S[,i] = c(min(X[,i]),max(X[,i])) } }
  if(n<par){ return(list(S,mean(y),n)) }
  if(length(y)<3){ return(list(S,mean(y),n)) }
  B = Bifurcate(X,y,S)
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

sum((Data[[2]]-PredictKrural3(B,Data[[1]]))^2)
B
Data[[1]]

B = Krural3(Data[[1]],Data[[2]],0.1)

B = Krural3(Xiris,Iris[,7],1)



plot(Xiris[,c(1,2)],col=1+Iris[,6])
