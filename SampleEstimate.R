SmpEstm <- function(Y,VZJ,theta,delta){
  # we have initial sample VZJ=list(V,Z,J,Jnum,Jsize,Jepoch) and theta
  # we are going to sample the posterior sample of {Jnum,Jsize,Jepoch}
  # and estimate the parameters
  N=length(Y)
  lambda=theta[1];mu=theta[2];beta=theta[3];a=theta[4];b=theta[5]
  if (is.na(lambda)) print(lambda)
  # First, construct a Markov Chain to Update {V_0,Z_0}
  JJJ=VZJ[4:6]
  V0=VZJ[[1]][1]
  for(i in 1:100){
    V0new=rgamma(1,a,b)
    accptP=V0accptProb(V0,V0new,Y,JJJ,theta,delta)
    if (runif(1)<accptP) V0=V0new
  }
  V=rep(0,N+1);Z=rep(0,N+1);V[1]=V0
  VZ=VZupdate(0,V,Z,JJJ,lambda,delta)
  V=VZ[[1]];Z=VZ[[2]]

  # Then, construct a Markov Chain to Update {Jnum[n],Jsize[[n]],Jepoch[[n]]}
  for(n in 1:N){
    for(i in 1:100){
      # propose
      Jnum_n=rpois(1,a*lambda*delta)
      if(Jnum_n>0){
        Jsize_n=rexp(Jnum_n,b);Jepoch_n=runif(Jnum_n,0,delta)
      } else {
        Jsize_n=0;Jepoch_n=0
      }
      JJJ_nnew=list(Jnum_n,Jsize_n,Jepoch_n)
      # update
      accptP=VnaccptProb(n,JJJ_nnew,V,Z,Y,JJJ,theta,delta)
      if (runif(1)<accptP) {
        JJJ[[1]][n]=Jnum_n
        JJJ[[2]][[n]]=Jsize_n
        JJJ[[3]][[n]]=Jepoch_n
      }
    }
    # update V,Z
    VZ=VZupdate(n,V,Z,JJJ,lambda,delta)
    V=VZ[[1]];Z=VZ[[2]]
  }

  # update the parameters lambda,mu,beta,a,b
  if (is.na(lambda)) print("lambda is na in SmpEstm")
  theta=ParamUpd(Y,JJJ,V,Z,lambda)
  return(theta)
}
