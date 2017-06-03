VZupdate <- function(n,V,Z,JJJ,lambda,delta){
  Jnum=JJJ[[1]];Jsize=JJJ[[2]];Jepoch=JJJ[[3]]
  N=length(V)-1
  # n=0, V[1],Z[1] need not to be changed
  if (n==0) n=n+1
  # n>0, it goes as following
  for(i in n:N){
    V[i+1]=exp(-lambda*delta)*V[i]
    Z[i+1]=Z[i]
    numJ=Jnum[i]
    if (numJ>0) {
      Z[i+1]=Z[i+1]+sum(Jsize[[i]])
      V[i+1]=V[i+1]+drop(exp((Jepoch[[i]]-delta)*lambda)%*%Jsize[[i]])
    }
  }
  VZ=list(V,Z)
  return(VZ)
}

LR <- function(Z,V,Znew,Vnew,Y,lambda,mu,beta){
  # compute the Likelihood Ratio
  diffZV=diff(Z)-diff(V)
  diffZVnew=diff(Znew)-diff(Vnew)
  Q=diffZV/lambda;MU=mu*delta+beta*Q
  Qnew=diffZVnew/lambda;MUnew=mu*delta+beta*Qnew
  part=-(Y-MUnew)^2/(2*Qnew)+(Y-MU)^2/(2*Q)
  logLR=(1/2)*sum(log(diffZV)-log(diffZVnew))+sum(part)

  return(exp(logLR))
}

V0accptProb<-function(V0,V0new,Y,JJJ,theta,delta){
  lambda=theta[1];mu=theta[2];beta=theta[3];a=theta[4];b=theta[5]
  # Jnum=VZJ[[4]];Jsize=VZJ[[5]];Jepoch=VZJ[[6]]
  N=length(Y)
  V=rep(0,N+1);V[1]=V0
  Vnew=rep(0,N+1);Vnew[1]=V0new;Z=rep(0,N+1)
  # Z is the same for V, Vnew
  VZ=VZupdate(0,V,Z,JJJ,lambda,delta)
  VZnew=VZupdate(0,Vnew,Z,JJJ,lambda,delta)
  V=VZ[[1]];Vnew=VZnew[[1]];Z=VZ[[2]];Znew=VZnew[[2]]
  # compute the likelihood ratio, lr
  lr=LR(Z,V,Znew,Vnew,Y,lambda,mu,beta)
  accptProb=min(lr,1)
  return(accptProb)
}

VnaccptProb<-function(n,JJJ_nnew,V,Z,Y,JJJ,theta,delta){
  # JJJ=list(Jnum,Jsize,Jepoch)
  # JJJ_nnew=list(Jnum_n,Jsize_n,Jepoch_n)
  N=length(Y)
  lambda=theta[1];mu=theta[2];beta=theta[3];a=theta[4];b=theta[5]

  JJJnew=JJJ
  JJJnew[[1]][n]=JJJ_nnew[[1]] # Jnum is a vector
  JJJnew[[2]][[n]]=JJJ_nnew[[2]] # Jsize is a list
  JJJnew[[3]][[n]]=JJJ_nnew[[3]] # Jepoch is a list
  VZnew=VZupdate(n,V,Z,JJJnew,lambda,delta)
  Vnew=VZnew[[1]];Znew=VZnew[[2]]

  ind=n:(N+1)
  lr = LR(Z[ind],V[ind],Znew[ind],Vnew[ind],Y[n:N],lambda,mu,beta)
  accptProb=min(lr,1)

  return(accptProb)
}
