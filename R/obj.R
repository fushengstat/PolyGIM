#### objective function for cc
obj <- function(x,y,u,rhoe,beta0,V,lambda,nu,tau,theta,alphas,betas,cases){
  k = length(lambda)
  n = length(y)
  
  Vi = solve(V)
  
  pr = prob(x,u,rhoe,lambda,nu,tau,theta,alphas,betas,cases)
  
  beta = listtovec(betas)$lst
  

  obj = sum(log(pr)) - 1/2*t(beta-beta0) %*% Vi %*% (beta-beta0)
 
  for (j in 1:k){
    yj = (y==j)
    aj = tau[j]
    bj = theta[,j]
    Dtj = Delta(x,aj,bj)
    
    obj = obj + sum(yj*log(Dtj))
  }
  

  obj = c(obj)
  return(obj)
}  

