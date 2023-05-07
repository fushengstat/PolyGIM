vars_ext <- function(models){
  nmodel = length(models)
  vars = rep(NA,nmodel)
  for(i in 1:nmodel){
    modeli = models[[i]]
    vars[i] = modeli$info$var
  }
  vars = unique(vars)

  vars
}

cases_ext <- function(models){
  nmodel = length(models)
  cases = rep(NA,nmodel)
  for(i in 1:nmodel){
    modeli = models[[i]]
    cases[i] = modeli$subtype
  }

  cases
}


polygim_v <- function(form,int, models, ncase, nctrl, V){
  formi = as.formula(form)
  nsample = list(ncase = ncase, nctrl = nctrl)


  vars = vars_ext(models)
  form0 = "y~"
  forme = as.formula(paste(form0, paste(vars, collapse="+")))

  rhoe = diag(ncase)/diag(nctrl)

  init0 = init(models,rhoe,int)
  alphas = init0$alphas
  betas = init0$betas
  beta0s = init0$beta0s

  beta0 = listtovec(beta0s)$lst

  y = int$y
  x <- model.matrix(formi, data = int)
  u <- model.matrix(forme, data = int)


  alphai = listtovec(alphas)$lst
  betai = listtovec(betas)$lst
  d.a = length(alphai)
  d.b = length(betai)


  nu = rep(0,d.a+d.b)

  lambda = lambda_value(y)
  rhoi = rho_value(y)
  k = length(lambda)
  d = ncol(x)-1

  int1 = int
  int1$y = int$y+1

  fit = vglm(formi, multinomial(refLevel = 1), data = int1)
  fita = summary(fit)
  tem1 = coef(fit, matrix = TRUE)
  tau0 = c(tem1[1,])-log(rhoi)
  theta0 = matrix(tem1[-1,],ncol =length(tau0))
  colnames(theta0) = names(coef(fit))[-(1:length(tau0))]

  nmodel = nrow(ncase)


  epsilon = 1e-3
  max_iter = 100
  # display = FALSE
  display = TRUE

  cases = cases_ext(models)
  fit1 = mygim(x,y,u,cases,rhoe,nsample,beta0,V,lambda,nu,tau0,theta0,alphas,
               betas,max_iter,epsilon,display)

  k = length(lambda)
  i1 = k + length(nu) + k
  i2 =  i1 + d*k

  theta1 = fit1$theta
  V1 = fit1$Var[(i1+1):i2,(i1+1):i2]
  nV1 = sqrt(diag(V1))
  se1 = matrix(nV1,d,k)

  colnames(theta1) = colnames(theta0)
  colnames(se1) = colnames(theta0)

  list(theta = theta1,se = se1)
}

polygim_opt <- function(form, int, models, ncase, nctrl){
  formi = as.formula(form)
  nsample = list(ncase = ncase, nctrl = nctrl)


  vars = vars_ext(models)
  form0 = "y~"
  forme = as.formula(paste(form0, paste(vars, collapse="+")))

  rhoe = diag(ncase)/diag(nctrl)

  init0 = init(models,rhoe,int)
  alphas = init0$alphas
  betas = init0$betas
  beta0s = init0$beta0s

  beta0 = listtovec(beta0s)$lst

  y = int$y
  x <- model.matrix(formi, data = int)
  u <- model.matrix(forme, data = int)


  alphai = listtovec(alphas)$lst
  betai = listtovec(betas)$lst
  d.a = length(alphai)
  d.b = length(betai)


  nu = rep(0,d.a+d.b)

  lambda = lambda_value(y)
  rhoi = rho_value(y)
  k = length(lambda)
  d = ncol(x)-1

  int1 = int
  int1$y = int$y+1

  fit = vglm(formi, multinomial(refLevel = 1), data = int1)
  a = summary(fit)
  tem1 = coef(fit, matrix = TRUE)
  tau0 = c(tem1[1,])-log(rhoi)
  theta0 = matrix(tem1[-1,],ncol =length(tau0))
  colnames(theta0) = names(coef(fit))[-(1:length(tau0))]


  nmodel = nrow(ncase)


  epsilon = 1e-3
  max_iter = 100
  # display = FALSE
  display = TRUE

  cases = cases_ext(models)
  fit1 = gim_iter(x,y,u,cases,rhoe,nsample,beta0,lambda,nu,tau0,theta0,alphas,
                  betas,max_iter,epsilon,display)
  k = length(lambda)
  i1 = k + length(nu) + k
  i2 =  i1 + d*k

  theta1 = fit1$theta
  V1 = fit1$Var[(i1+1):i2,(i1+1):i2]
  nV1 = sqrt(diag(V1))
  se1 = matrix(nV1,d,k)

  colnames(theta1) = colnames(theta0)
  colnames(se1) = colnames(theta0)
  list(theta = theta1,se = se1)
}


