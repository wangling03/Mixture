PFAdecomp <- function (hat_beta, Sigma, smallestf,e, stdized) 
{

  hat_beta <- as.vector(hat_beta)
  Sigma <- as.matrix(Sigma)
  
  P <- length(hat_beta)
  
  # standardize the data
  
  SD <- sqrt(diag(Sigma))
  
  switch(stdized, 
         "Y"={  hat_beta <- hat_beta/SD
                Sigma <- diag(1/SD) %*% Sigma %*% diag(1/SD)
         },
         "N"={ hat_beta <- hat_beta
               Sigma <- Sigma 
         }
  )
  
  #Spectral Decomposition of the covariance matrix
  
  # pca 
  Kmax=P

  Sigma_spec=eigen(Sigma, symmetric=TRUE, only.values = FALSE, EISPACK = FALSE)
  lambda=Sigma_spec$values
  eigvec=Sigma_spec$vectors
  
  # determine the factor loadings
  K=0
  
  if (K==0) {
    K = 1
    while( (K < Kmax) & (sqrt(sum(lambda[(K+1):length(lambda)]^2)) >= e*sum(lambda))) 
      K = K + 1 
  }    
  
  print(K)
  
  sqrt_lambda <- as.matrix(sqrt(lambda[1:K]))
  b <- as.matrix(eigvec[,1:K])
  for (i in 1:K)  {
    b[,i] <- b[,i] * sqrt_lambda[i]  # factor loadings
  }
  
  # regression
  # estimate the factors, with 10% largest |Z| eliminated.
  
  o = order(abs(hat_beta))
  Zperm = hat_beta[o]
  Lperm = as.matrix(b[o,])
  Z.reduce = Zperm[1:(round(P*smallestf))]
  L.reduce = as.matrix(Lperm[1:(round(P*smallestf)),]) 
  
  W.hat=L1linreg(L.reduce, Z.reduce, p = 1, tol = 1e-07, maxiter = 5000)$x
  
  #W.hat <- rq(Z.reduce~L.reduce-1, 0.5)$coef

  #hat_beta_rlm<-rlm(Z.reduce~L.reduce-1, maxit=500)
  #W.hat<-numeric(0)
  #W.hat<-hat_beta_rlm$coefficients
  
  
  #a1<-lmrob.control(k.max = 1000)
  #hat_beta_lmrobfit <- lmrob(Z.reduce~L.reduce-1,control=a1)
  #W.hat<-numeric(0)
  #W.hat<-hat_beta_lmrobfit$coefficients
  
  bW.est <- b%*%(W.hat)
  
  z=hat_beta-bW.est
  
  #Calculating variance for W
  var_K=numeric(0)
  for (j in 1:P)  
  {
    var_K[j]=Sigma[j,j]-sum(b[j,]^2)
  }      
  
  output<-list(z,var_K,SD)

  return(output)
}
