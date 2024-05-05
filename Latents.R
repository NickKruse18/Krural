









Equation = function(x){
  return(x[1]*x[2])
}

MatrixEquation = function(n){
  X = matrix(0,n,2)
  y = numeric(n)
  for(i in 1:n){
    X[i,] = runif(2,0,10)
    y[i] = Equation(X[i,])
  }
  return(list(X,y))
}


Data = MatrixEquation(10000)

summary(lm(Data[[2]] ~ Data[[1]]))

B = KruralNickwork(Data[[1]],Data[[2]],0.00001)

sum((PredictKruralNickwork(B,Data[[1]]) - Data[[2]])^2)

cbind(Data[[1]],Data[[1]][,1]*Data[[1]][,2],PredictKruralNickwork(B,Data[[1]]))





