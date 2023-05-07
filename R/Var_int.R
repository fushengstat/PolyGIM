## Variance of internal study--theta
var_int <- function(x,y,u,rhoe,nsample,V,lambda,nu,tau,theta,alphas,betas,cases){
  n = nrow(x)
  k = length(lambda)
  H = hess(x,y,u,rhoe,V,lambda,nu,tau,theta,alphas,betas,cases)
  td = nrow(H)
  
  alpha = listtovec(alphas)$lst
  beta = listtovec(betas)$lst


  ## first term
  nH = matrix(0,td,td)
  H0 = -H/n
  i2 = length(lambda)+length(nu)
  i3 = td - length(beta)
  nH[1:i2,1:i2] = -H0[1:i2,1:i2]
  nH[(i2+1):td,(i2+1):td] = H0[(i2+1):td,(i2+1):td]
  jbeta = nH[(i3+1):td,(i3+1):td]
  
  
  ## second term 
  J1 = H0[,1:k]
  g1 = matrix(lambda,k,1)
  Gm = diag(lambda) - g1%*%t(g1)
  
  ## third term 
  d2 = length(alpha)
  dt = d2 + length(beta)
  H1 = Sigma0_ext(x,u,rhoe,nsample,lambda,nu,tau,theta,alphas,betas,cases)
  Sig0 = H1[(d2+1):dt,(d2+1):dt]
  gam = 1/n
  
  Q = matrix(0,td,td)
  Vi = solve(V)
  Q[(i3+1):td,(i3+1):td] = gam*Vi%*%Sig0%*%Vi - jbeta
  
  
  I0 = nH - J1 %*% Gm %*% t(J1) + Q
  
  nS = 1/n* solve(H0) %*% I0 %*% solve(H0) 
  
  return(nS)
}

