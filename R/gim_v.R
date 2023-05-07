#### polygim procedure with fixed V
mygim <- function(x,y,u,cases,rho,nsample,beta0,V,lambda,nu,tau,theta,alphas,betas,
                  max_iter,epsilon,display){

  alpha = listtovec(alphas)$lst
  id.alp = listtovec(alphas)$leng
  d2 = sum(id.alp)
  da.all = getindex(id.alp)
  beta = listtovec(betas)$lst
  id.bet = listtovec(betas)$leng
  d3 = sum(id.bet)
  db.all = getindex(id.bet)

  nm.alp = rownames(alpha)
  nm.bet = rownames(beta)


  d1 = ncol(x)-1
  d2 = length(alpha)
  d3 = length(beta)
  k = length(lambda)


  ## some index sets
  i1 = k
  i2 = i1 + d2 + d3
  i3 = i2 + k
  i4 = i3 + k*d1
  i5 = i4 + d2
  td = i5 + d3


  err = Inf
  iter = 1

  tobj = rep(0,max_iter)
  tobj[1] = obj(x,y,u,rho,beta0,V,lambda,nu,tau,theta,alphas,betas,cases)

  while ((err>=epsilon)&(iter<=max_iter)){
    # ## old values of parameters
    ttheta = c(theta)
    omu = c(lambda,nu,tau,ttheta,alpha,beta)

    ## new parameters updating via newton-raphson algorithm
    gd = grad(x,y,u,rho,beta0,V,lambda,nu,tau,theta,alphas,betas,cases)
    H = hess(x,y,u,rho,V,lambda,nu,tau,theta,alphas,betas,cases)
    nmu = omu - solve(H)%*%gd

    ### new parameters
    lambda = nmu[1:i1]
    nu = nmu[(i1+1):i2]
    tau = nmu[(i2+1):i3]
    theta = matrix(nmu[(i3+1):i4],nrow = d1,ncol = k)
    alpha = nmu[(i4+1):i5]
    beta = nmu[(i5+1):td]

    names(alpha) = nm.alp
    names(beta) = nm.bet

    err = norm(nmu-omu, type="2")
    iter = iter+1
    alphas = vectolist(alpha,da.all)
    betas = vectolist(beta,db.all)

    tobj[iter] = obj(x,y,u,rho,beta0,V,lambda,nu,tau,theta,alphas,betas,cases)

    if (display){
      print(paste("Iteration ", iter, ", err: ", err, sep =""))
    }

  }

  ### summary
  tobj = tobj[1:iter]

  V1 = var_int(x,y,u,rho,nsample,V,lambda,nu,tau,theta,alphas,betas,cases)

  list(lambda = lambda, nu = nu,tau = tau, theta = theta,
       alphas = alphas,betas = betas,Var = V1, obj = tobj)

}


