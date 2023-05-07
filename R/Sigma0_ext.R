## Variance of Sigma_0 iteration 
Sigma0_ext <- function(x,u,rhoe,nsample,lambda,nu,tau,theta,alphas,betas,cases){
  n = nrow(x)
  pr = prob(x,u,rhoe,lambda,nu,tau,theta,alphas,betas,cases)
  pr = pr/n
  
  ncase = nsample$ncase
  nctrl = nsample$nctrl
  
  K = length(betas)
  
  alpha = listtovec(alphas)$lst
  id.alp = listtovec(alphas)$leng
  beta = listtovec(betas)$lst
  id.bet = listtovec(betas)$leng
  de.all = getindex(id.alp+id.bet)
  index.ext  = mergeid(id.alp,id.bet)
  
  ida = index.ext$id1
  idb = index.ext$id2
  
  
  d0 = length(nu)
  nV = matrix(0,d0,d0)
  
  A = matrix(0,d0,d0)
  B = matrix(0,d0,d0)
  w0 = matrix(0,d0,1)
  
  for (s1 in 1:K){ 
    alphaj = alphas[[s1]]
    betaj = betas[[s1]]
    c1 = cases[s1]
    
    names = colnames(u)
    nm1 = c(names(alphaj),names(betaj))
    u1 = u[,names %in% nm1]
    
    rho1 = rhoe[s1]
    p01_1s = phis(u1,rho1,alphaj,betaj)
    p0_1 = p01_1s$p0
    p1_1 = p01_1s$p1
    
    
    a1 = tau[c1]
    b1 = theta[,c1]
    Dt1 = Delta(x,a1,b1)
    dt1 = delta(u1,alphaj,betaj)
    
    is1 = index.ext$id3[[s1]]
    n1.ss = ncase[s1,s1]
    
    ## hessian matrix
    t1 = - n1.ss*pr*dt1*(1+rho1*Dt1)/((1+rho1*dt1)^2)
    A[is1,is1] = t(c(t1)*u1)%*% u1
    w0[is1] = 1/rho1*t(p0_1)%*%pr
    
    for(s2 in 1:K){
      alpha2 = alphas[[s2]]
      beta2 = betas[[s2]]
      c2 = cases[s2]
      n0.ij = nctrl[s1,s2]
      n1.ij = ncase[s1,s2]
      
      
      names = colnames(u)
      nm2 = c(names(alpha2),names(beta2))
      u2 = u[,names %in% nm2]
      
      rho2 = rhoe[s2]
      p01_2s = phis(u2,rho2,alpha2,beta2)
      p0_2 = p01_2s$p0
      p1_2 = p01_2s$p1
      
      a2 = tau[c2]
      b2 = theta[,c2]
      Dt2 = Delta(x,a2,b2)
      
      is2 = index.ext$id3[[s2]]
      
      ## covariance matrix
      B[is1,is2] = n0.ij*t(c(pr)*p0_1)%*%p0_2+n1.ij*t(c(Dt2*pr)*p1_1)%*%p1_2
    }
  }
  
  for(s in 1:K){
    is = index.ext$id3[[s]]
    rhos = rhoe[s]
    ws = w0[is]
    for(t in 1:K){
      it = index.ext$id3[[t]]
      n0.ij = nctrl[s,t]
      n1.ij = ncase[s,t]
      rhot = rhoe[t]
      
      wt = w0[it]
      B[is,it] = B[is,it]- (n1.ij+n0.ij*rhot*rhos)*ws%*%(t(wt))
    }
  }
  
  Sig0 = solve(A)%*%B%*%solve(A)
  
  nid = c(ida,idb)
  Sig = Sig0[nid,nid]

  return(Sig)
}

