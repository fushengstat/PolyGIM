#### hessian matrix function
hess <- function(x,y,u,rhoe,V,lambda,nu,tau,theta,alphas,betas,cases){
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
  id1 = index.ext$id1
  id2 = index.ext$id2

  ## some index sets
  i1 = k
  i2 = i1 + d2 + d3
  i3 = i2 + k
  i4 = i3 + k*d1 
  i5 = i4 + d2
  td = i5 + d3
  
  ### hessian
  h0 = matrix(0,td,td)
  h1 = matrix(0,td,td)
  

  
  ga = g_all(x,u,rhoe,tau,theta,alphas,betas,cases)
  g_nu = c(pr)*ga
  
  g_lam = matrix(0,n,k)
  g_tau = matrix(0,n,k)
  g_theta = list()
  for(j in 1:k){
    aj = tau[j]
    bj = theta[,j]
    Dtj = Delta(x,aj,bj)
    
    g_lam[,j] = (Dtj-1)*pr
    g_tau[,j] = lambda[j]*pr*Dtj
    g_theta[[j]] = lambda[j]*c(pr*Dtj)*x[,-1]
  }
  
  
  g_ext = matrix(0,n,d2+d3)
  for (s in 1:K){
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    uj = u[,names %in% nmj]
    
    rhoj = rhoe[s]
    p1 = phi1(uj,rhoj,alphaj,betaj)

    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    dt = delta(uj,alphaj,betaj)
    
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    
    t0 =  rhoj*(p1%*%nus)
    g_tau[,cj] = g_tau[,cj] + t0*pr*Dt 
    g_theta[[cj]] = g_theta[[cj]] + c(t0*pr*Dt)*x[,-1]
    
    
    t2 = uj%*%nus
    tm = -rhoj*pr*t2*(1+rhoj*Dt)/((1+rhoj*dt)^2)*dt
    g_ext[,nuid] = c(tm)*uj
  }
  
  g0 = matrix(0,n,td)
  g0[,1:i1] = g_lam
  g0[,(1+i1):i2] = g_nu
  g0[,(1+i2):i3] = g_tau
  g0[,(i4+1):i5] = g_ext[,id1]
  g0[,(i5+1):td] = g_ext[,id2]
  for (j in 1:k){
    n.id1 = i3+(j-1)*d1+1
    n.id2 = i3+j*d1
    g0[,n.id1:n.id2] = g_theta[[j]]
  }
  
  h0  = t(g0)%*%g0 
  
  ## block 1: partial derivative: lambda - vs others (line 1:k---i1)
  for (s in 1:k){ 
    aj = tau[s]
    bj = theta[,s]
    Dtj = Delta(x,aj,bj)
    
    h1[s,i2+s] = - t(pr)%*%Dtj
    h1[i2+s,s] = - t(pr)%*%Dtj
    
    j1 = i3+(s-1)*d1+1
    j2 = i3+s*d1
    h1[s,j1:j2] = - t(pr*Dtj)%*%x[,-1]
    h1[j1:j2,s] = - t(pr*Dtj)%*%x[,-1]
  }
  
  ## block 2: partial derivative: nu- vs others (line i1+1---i2)
  for (s in 1:K){ 
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    uj = u[,names %in% nmj]
    
    rhoj = rhoe[s]
    p1 = phi1(uj,rhoj,alphaj,betaj)

    
    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    dt = delta(uj,alphaj,betaj)
    
    tm0 = rhoj*c(pr*Dt)*p1
    h1[(i1+nuid),i2+cj] = - colSums(tm0)
    h1[i2+cj,(i1+nuid)] = - colSums(tm0)
    
    tm1 = rhoj*c(Dt)*p1
    h1[(i1+nuid),i2+cj] = - t(pr) %*% tm1
    h1[i2+cj,(i1+nuid)] = - colSums(tm0)
    
    
    
    j1 = i3+(cj-1)*d1+1
    j2 = i3+cj*d1
    tm1 = t(tm0)%*%x[,-1]
    h1[(i1+nuid),j1:j2] = - tm1
    h1[j1:j2,(i1+nuid)] = - t(tm1)
    
    
    sj = index.ext$id3[[s]]
    tem = -pr*rhoj*(1+rhoj*Dt)/((1+rhoj*dt)^2)*dt
    ntem = c(tem)*uj
    h1[i1+sj,i4+sj] = -t(ntem) %*% uj
    h1[i4+sj,i1+sj] = -t(ntem) %*% uj
  }
  
  ## block 3: partial derivative: tau - vs others (line (i2+1)---i3)
  for (s in 1:k){ 
    aj = tau[s]
    bj = theta[,s]
    Dtj = Delta(x,aj,bj)
    
    i0 = i2+s
    h1[i0,i0] = -lambda[s]*t(Dtj)%*%pr
    
    j1 = i3+(s-1)*d1+1
    j2 = i3+s*d1
    h1[i0,j1:j2] = - lambda[s]*t(pr*Dtj)%*%x[,-1]
    h1[j1:j2,i0] = - lambda[s]*t(pr*Dtj)%*%x[,-1]
  }
  
  for(s in 1:K){ 
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    uj = u[,names %in% nmj]
    
    rhoj = rhoe[s]
    p1 = phi1(uj,rhoj,alphaj,betaj)

    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    dt = delta(uj,alphaj,betaj)
    
    i0 = i2+cj
    j1 = i3+(cj-1)*d1+1
    j2 = i3+cj*d1
    t1 = p1 %*% nus
    h1[i0,i0] = h1[i0,i0] - rhoj*t(Dtj*pr)%*%t1
    h1[i0,j1:j2] = h1[i0,j1:j2] - rhoj*t(Dtj*pr*t1)%*%x[,-1]
    h1[j1:j2,i0] = h1[j1:j2,i0] - rhoj*t(Dtj*pr*t1)%*%x[,-1]
    
    t2 = uj%*%nus
    t3 = rhoj*Dt
    t4 = -rhoj*pr*t2*t3/((1+rhoj*dt)^2)*dt
    tem1 = t(t4)%*%uj
    sj = index.ext$id3[[s]]
    h1[i0,(i4+sj)] = h1[i0,(i4+sj)] - tem1
    h1[(i4+sj),i0] = h1[(i4+sj),i0] - t(tem1)
  }
  
  
  ## block 4: partial derivative: theta - vs others (line (i3+1)---i4)
  for (s in 1:k){ 
    aj = tau[s]
    bj = theta[,s]
    Dtj = Delta(x,aj,bj)
    
    j1 = i3+(s-1)*d1+1
    j2 = i3+s*d1
    h1[j1:j2,j1:j2] = - lambda[s]*t(c(pr*Dtj)*x[,-1])%*%x[,-1]
  }
  
  for (s in 1:K){ 
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    uj = u[,names %in% nmj]
    
    
    rhoj = rhoe[s]
    p1 = phi1(uj,rhoj,alphaj,betaj)

    
    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    dt = delta(uj,alphaj,betaj)
    
    
    
    j1 = i3+(cj-1)*d1+1
    j2 = i3+cj*d1
    t1 = p1%*%nus
    t2 = c(t1*pr*Dtj)
    h1[j1:j2,j1:j2] = h1[j1:j2,j1:j2] - rhoj*t(t2*x[,-1])%*%x[,-1]
    
    t3 = uj%*%nus
    t4 = rhoj*Dt
    t5 = rhoj*pr*t3*t4/((1+rhoj*dt)^2)*dt
    tem2 = -t(c(t5)*x[,-1])%*%uj
    sj = index.ext$id3[[s]]
    h1[j1:j2,(i4+sj)] = h1[j1:j2,(i4+sj)] - tem2
    h1[(i4+sj),j1:j2] = h1[(i4+sj),j1:j2] - t(tem2)      
  }
  
  
  ## block 5: partial derivative: ext - vs others (line (i4+1)---td)
  for(s in 1:K){
    alphaj = alphas[[s]]
    betaj = betas[[s]]
    cj = cases[s]
    
    names = colnames(u)
    nmj = c(names(alphaj),names(betaj))
    uj = u[,names %in% nmj]
    
    rhoj = rhoe[s]
    p1 = phi1(uj,rhoj,alphaj,betaj)
    a1 = tau[cj]
    b1 = theta[,cj]
    Dt = Delta(x,a1,b1)
    nuid = de.all$sid[s]:de.all$eid[s]
    nus = nu[nuid]
    dt = delta(uj,alphaj,betaj)
    
    t1 = uj%*%nus
    t2 = -rhoj*pr*(1+rhoj*Dt)*dt*(1-rhoj*dt)/((1+rhoj*dt)^3)
    tem3 = t(c(t1*t2)*uj)%*% uj
    
    is = i4+index.ext$id3[[s]]
    h1[is,is] = h1[is,is]-tem3
  }

  ida = i4 + index.ext$id1
  idb = i4 + index.ext$id2
  h1[idb,idb] = h1[idb,idb]-Vi
  
  
  id.new = c(1:i4,ida,idb)
  h1 = h1[id.new,id.new]
  
  H = h0+h1
  
  return(H)
}


