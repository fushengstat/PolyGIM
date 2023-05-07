#### gradient vector function
grad <- function(x,y,u,rhoe,beta0,V,lambda,nu,tau,theta,alphas,betas,cases){
  k = length(lambda)
  Vi = solve(V)
  n = length(y)
  K = length(cases)
  
  pr = prob(x,u,rhoe,lambda,nu,tau,theta,alphas,betas,cases)
  

  d1 = ncol(x)-1
  dnu = length(nu)
  ## dnu = d2+d3
  alpha = listtovec(alphas)$lst
  id.alp = listtovec(alphas)$leng
  d2 = sum(id.alp)
  da.all = getindex(id.alp)
  beta = listtovec(betas)$lst
  id.bet = listtovec(betas)$leng
  d3 = sum(id.bet)
  db.all = getindex(id.bet)
  
  de.all = getindex(id.alp+id.bet)
  index.ext  = mergeid(id.alp,id.bet)

  
  g = g_all(x,u,rhoe,tau,theta,alphas,betas,cases)
  grad_nu = - t(g) %*% pr
  
  grad_lam = rep(0,k)
  grad_tau= rep(0,k)
  grad_theta = matrix(0,d1,k)
  for(j in 1:k){
    aj = tau[j]
    bj = theta[,j]
    Dtj = Delta(x,aj,bj)
    idj = (y==j)
    
    lamj = lambda[j]
    grad_lam[j] = -t(Dtj-1) %*% pr
    
    grad_tau[j] = sum(idj) - lamj*(t(Dtj)%*% pr) 
    grad_theta[,j] = t(idj)%*%x[,-1] - lamj*(t(Dtj*pr)%*%x[,-1])
  }
  
  grad_ext = rep(0,d2+d3)
  for (s in 1:K){
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    Uj = u[,names %in% nmj]
    
    rhoj = rhoe[s]
    p01 = phis(Uj,rhoj,alphaj,betaj)
    p1 = p01$p1
    
    dt = delta(Uj,alphaj,betaj)
    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    
    
    t0 =  rhoj*(p1%*%nus)
    grad_tau[cj] = grad_tau[cj] - sum(t0*pr*Dt)
    grad_theta[,cj] = grad_theta[,cj] - t(x[,-1]) %*% (t0*pr*Dt)
    
    t1 = Uj %*% nus
    tem = -t1*rhoj*dt*(1+rhoj*Dt)/((1+rhoj*dt)^2)
    grad_ext[nuid] =  - t(Uj)%*% (pr*tem)
  }

  id1 = index.ext$id1
  id2 = index.ext$id2
  
  grad_alpha = grad_ext[id1]
  grad_beta = grad_ext[id2]- Vi%*%(beta-beta0)
  
  
  gtheta = c(grad_theta)
  gradient = c(grad_lam,grad_nu,grad_tau,gtheta,grad_alpha,grad_beta)
  
  ngrad = matrix(gradient,ncol=1)
  
  return(ngrad)
}

