lambda_value <- function(y) {
  uy <- unique(y)
  k <- length(uy) - 1
  n <- length(y)
  
  lambda <- integer(k)
  counts <- tabulate(y)
  lambda[1:k] <- counts[1:k] / n
  
  return(lambda)
}

rho_value <- function(y) {
  uy <- unique(y)
  k <- length(uy) - 1
  
  counts <- table(y)
  rho <- sapply(2:(k+1), function(j) counts[j] / counts[1])
  
  return(rho)
}



##### ----linear---external model 
delta <- function(U,a,b) {
  obj <- exp(U %*% c(a,b))
  return(obj)
}


## phi0=
phi0 <- function(U,rhoe,alpha,beta){
  dt = delta(U,alpha,beta)
  a1 = -(rhoe*dt)/(1+rhoe*dt)
  p0 = c(a1) * U
  
  return(p0)
}


## phi1=
phi1 <- function(U,rhoe,alpha,beta){
  dt = delta(U,alpha,beta)
  a1 = 1/(1+rhoe*dt)
  p1 = c(a1) * U
  
  return(p1)
}

phis <- function(U, rhoe, alpha, beta) {
  dt <- delta(U, alpha, beta)
  a1 <-  1/ (1 + rhoe * dt)
  a0 <- -(rhoe * dt) *a1
  p0 <- sweep(U, 1, a0, "*")
  
  p1 <- sweep(U, 1, a1, "*")
  
  
  return(list(p0 = p0, p1 = p1))
}




##### ----linear-- internal model
Delta <- function(X,a,b){
  obj = exp(X %*% c(a,b))
  return(obj)
}






## external data -- estimation equation for control vs case_j
gj.func <- function(X,Uj,rhoj,tauj,thetaj,alphaj,betaj){
  p01 = phis(Uj,rhoj,alphaj,betaj)
  p0 = p01$p0
  p1 = p01$p1
  
  
  Dtj = Delta(X,tauj,thetaj)
  Dtj = as.vector(Dtj)
  
  gj = p0 + rhoj*(Dtj*p1)
  
  return(gj)
}



## all vectors for estimation equation
g_all <- function(X,U,rhoe,tau,theta,alphas,betas,cases){
  K = length(cases)
  
  gall = NULL
  
  for(i in 1:K){
    cj = cases[i]
    alphaj = alphas[[i]]
    betaj = betas[[i]]
    
    rhoj = rhoe[i]
    tauj = tau[cj]
    thetaj = theta[,cj]
    
    
    ## it does not contain the intercept
    names = colnames(U)
    nmj = c(names(alphaj),names(betaj))
    Uj = U[,names %in% nmj]
    
    gj = gj.func(X,Uj,rhoj,tauj,thetaj,alphaj,betaj)
    
    gall = cbind(gall,gj)
  }
  
  return(gall)
}


## prob for case_j vs control
prob <- function(X,U,rhoe,lambda,nu,tau,theta,alpha,beta,cases){
  n = nrow(X)
  k = length(lambda)
  

  np = rep(0,n)
  for (j in 1:k){
    aj = tau[j]
    bj = theta[,j]
    Dtj = Delta(X,aj,bj)
    
    np = np + lambda[j]*(Dtj-1)
  }
  
  g = g_all(X,U,rhoe,tau,theta,alpha,beta,cases)
  
  p  = 1/(1+np+g%*%nu)
  
  return(p)
}



## listtovec
listtovec <- function(lists) {
  K <- length(lists)
  lst <- unlist(lists, use.names = FALSE)
  id <- lengths(lists)
  nms <- names(unlist(lists))
  lst <- matrix(lst, ncol = 1, dimnames = list(nms, NULL))
  list(lst = lst, leng = id)
}




vectolist <- function(vecs,lens){
  K = length(lens$eid)
  lst = list()
  for(i in 1:K){
    id = lens$sid[i]:lens$eid[i]
    lst[[i]] = vecs[id]
    
    names(lst[[i]]) = names(vecs)[id]
  }
  lst
}


getindex <- function(lens){
  K = length(lens)
  sid = c(1,1+cumsum(lens))[1:K]
  eid = cumsum(lens)
  
  list(sid = sid,eid = eid)
}




mergeid <- function(len1,len2){
  K = length(len1)
  tid = len1+len2
  id.all = getindex(tid)
  
  id1 <-NULL
  id2 <-NULL
  id3 <- list()

  
  for(j in 1:K){
    s1 = id.all$sid[j]
    e1 = s1+len1[j]-1
    s2 = e1+1
    e2 = id.all$eid[j]
    
    id1 = c(id1,s1:e1)
    id2 = c(id2,s2:e2)
    id3[[j]] = c(s1:e2)
  }

  list(id1 = id1, id2 = id2,id3 = id3)
  
}




init <- function(models,rhoe,int){
  K = length(models)
  
  alphas = list()
  betas = list()
  beta0s = list()
  Se0s = list()
  for(i in 1:K){
    formi = models[[i]][[1]]
    casei = models[[i]][[4]]
    
    valp =  models[[i]][[2]]
    vbet =  models[[i]][[3]]$var
    beta0s[[i]] = models[[i]][[3]]$bet
    names(beta0s[[i]]) = models[[i]][[3]]$var
    
    Se0s[[i]] = models[[i]][[3]]$se
    names(Se0s[[i]]) = models[[i]][[3]]$var
    
    
    
    f <- as.formula(formi)
    mat <- model.matrix(f, data = int)
    
    group = c(0,casei)
    
    ni = sum(int$y==casei)
    n0 = sum(int$y==0)
    
    
    id = which(int$y %in% group)
    int0 = int[id,]
    int0$y = as.numeric(int0$y!=0)
    fit <- glm(f, family = 'binomial', data = int0)
    
    alphas[[i]] = coef(fit)[valp]
    alphas[[i]][1] =  alphas[[i]][1] - log(ni/n0) + log(rhoe[i])
    betas[[i]] = coef(fit)[vbet]
    
  }
  
  list(alphas=alphas,betas = betas, beta0s = beta0s, Se0s = Se0s)
}

