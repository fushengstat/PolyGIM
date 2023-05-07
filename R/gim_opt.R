# iterative procedure of ploygim
gim_iter <- function(x,y,u,cases,rho,nsample,beta0,lambda,nu,tau,theta,alphas,betas,
                     max_iter,epsilon,display){
  alpha = listtovec(alphas)$lst
  id.alp = listtovec(alphas)$leng
  d2 = sum(id.alp)
  da.all = getindex(id.alp)
  beta = listtovec(betas)$lst
  id.bet = listtovec(betas)$leng
  d3 = sum(id.bet)
  db.all = getindex(id.bet)


  d1 = ncol(x)
  d2 = length(alpha)
  d3 = length(beta)

  ## some index sets
  dt = d2 + d3

  err = Inf
  iter = 1

  tobj = rep(0,max_iter)

  while ((err>=epsilon)&(iter<=max_iter)){
    ttheta = c(theta)
    omu = c(lambda,nu,tau,ttheta,alpha,beta)

    ## initial variance of beta
    H1 = Sigma0_ext(x,u,rho,nsample,lambda,nu,tau,theta,alphas,betas,cases)
    V0 = H1[(d2+1):dt,(d2+1):dt]

    if(iter == 1){
      tobj[1] =obj(x,y,u,rho,beta0,V0,lambda,nu,tau,theta,alphas,betas,cases)
    }
    ## main gim procedure
    res = mygim(x,y,u,cases,rho,nsample,beta0,V0,lambda,nu,tau,theta,
                alphas,betas,max_iter,epsilon,display = FALSE)


    ### new parameters
    lambda = res$lambda
    nu = res$nu
    tau = res$tau
    theta = res$theta
    alphas = res$alpha
    betas = res$beta
    ttheta = c(theta)

    alpha = listtovec(alphas)$lst
    beta = listtovec(betas)$lst


    nmu = c(lambda,nu,tau,ttheta,alpha,beta)

    err = norm(nmu-omu, type="2")
    iter = iter+1

    tobj[iter] =obj(x,y,u,rho,beta0,V0,lambda,nu,tau,theta,alphas,betas,cases)

    if (display){
      print(paste("Iteration ", iter, ", err: ", err, sep =""))
    }
  }


  ### summary of variance
  H1 = Sigma0_ext(x,u,rho,nsample,lambda,nu,tau,theta,alphas,betas,cases)
  V0 = H1[(d2+1):dt,(d2+1):dt]
  V1 = var_int(x,y,u,rho,nsample,V0,lambda,nu,tau,theta,alphas,betas,cases)


  tobj = tobj[1:iter]

  list(lambda = lambda,nu = nu,tau = tau, theta = theta,
       alphas = alphas,betas = betas, obj = tobj, Var = V1)

}





