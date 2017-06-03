SmpGammaOU <- function(a,b,lambda,delta,N){
  V <- rep(0,N); Z <- rep(0,N)
  V[1] <- rgamma(1,shape = a,rate = b)
  J <- rep(FALSE,N-1)
  Jnum <- rpois(N-1,a*lambda*delta); all <- sum(Jnum)
  Jepoch <- as.list(rep(0,N-1)) # Jepoch <- vector("list",N-1)
  Jsize <- as.list(rep(0,N-1)) # Jsize <- vector("list",N-1)

  T_all <- runif(all,0,delta)
  Jsize_all <- rexp(all,b)

  ind <- 0

  for(n in 2:N){
    V[n] = exp(-lambda*delta)*V[n-1]
    Z[n]=Z[n-1]
    i=Jnum[n-1]
    if (i != 0) {
      J[n-1]=TRUE
      T_n=T_all[(ind+1):(ind+i)]
      Jsize_n=Jsize_all[(ind+1):(ind+i)]
      Jsize[[n-1]]=Jsize_n
      Jepoch[[n-1]]=T_n
      Z[n]=Z[n]+sum(Jsize_n)
      V[n]=V[n]+drop(exp((T_n-delta)*lambda)%*%Jsize_n)
      # %*% results in a matrix, use drop()
      ind=ind+i
    }
  }
  result=list(V,Z,J,Jnum,Jsize,Jepoch)
  return(result)
}

SampleY <- function(theta,delta,N){
  lambda = theta[1];a=theta[4];b=theta[5]
  mu=theta[2];beta=theta[3]
  VZ <- SmpGammaOU(a,b,lambda,delta,N)
  V=VZ[[1]];Z=VZ[[2]]
  Q=(diff(Z)-diff(V))/lambda
  Y=rnorm(N-1,mu*delta+beta*Q,sqrt(Q))
  return(list(Y,VZ))
}
