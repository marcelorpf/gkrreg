sigma2est1 = function(y1, y2, frac = .5) {
  n = length(y1)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmp = (y1[idx1] - y2[idx2])^2
  mean(quantile(tmp[tmp != 0], probs = c(.9, .1)))
}

sigma2est2 = function(y1, y2) {
  n = length(y1)
  tmp = c()
  for (i in 1:n) {
    for (j in 1:n) {
      tmp = c(tmp, (y1[i] - y2[j])^2)
    }
  }
  tmp = tmp[tmp != 0]
  median(tmp)
}

###### KERNEL FUNCTIONS ###############################

gauss.kern = function(a, b, s)
{
  #  as.vector(exp(-(1/sqrt(s))*abs(a-b))) #Laplacian 
  as.vector(exp(-(1/s)*(a-b)^2)) #Gaussian
  #  as.vector(exp(-(1/s)*abs(a-b))) #Exponential
}

#### KRR METHODS - ACCORDING TO THE HYPER-PARAMETER ####

kernel.reg1 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sigma2est1(y, yhat)  #################
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}

kernel.reg2 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sigma2est2(y, yhat) #################
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}

kernel.reg6 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  s2 = sum((y-yhat)^2)/(n-ncol(x)) #################
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}
