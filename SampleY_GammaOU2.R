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
      V[n]=V[n]+exp((T_n-delta)*lambda)%*%Jsize_n
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
  return(Y)
}

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
      V[i+1]=V[i+1]+exp((Jepoch[[i]]-delta)*lambda)%*%Jsize[[i]]
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

lmbdGrad <- function(Y,JJJ,V,Q,mu,beta,lambda){
  Jnum=JJJ[[1]];Jsize=JJJ[[2]];Jepoch=JJJ[[3]]
  N=length(Jnum);delta=1
  # Q_n=((Z_n-Z_{n-1})-(V_n-V_{n-1}))/lambda
  # V_n=exp(-lambda*delta)+\sum_1^{N(\lambda\Delta)}exp((T_j-delta)
  # *lambda)*J_j
  Vlmbd=rep(0,N+1)
  for(i in 1:N){
    numJ=Jnum[i]
    Vlmbd[i+1]=exp(-lambda*delta)*(Vlmbd[i]-delta*V[i])
    if(numJ>0){
      epk=Jepoch[[i]]-delta
      Vlmbd[i+1]=Vlmbd[i+1]+(exp(epk*lambda)*epk)%*%Jsize[[i]]
    }
  }
  Qlmbd=-(diff(Vlmbd)+Q)/lambda
  devi=Y-(mu*delta+beta*Q)
  logGlmbd=((-1/2+beta*devi)/Q+(1/2)*(devi/Q)^2)%*%Qlmbd
  return(as.numeric(logGlmbd))
}

ParamUpd <- function(Y,JJJ,V,Z,lambda){
  Jnum=JJJ[[1]];Jsize=JJJ[[2]];Jepoch=JJJ[[3]]
  N=length(Y);delta=1
  # estimate b
  meanJsize=sum(sapply(Jsize,sum))/sum(Jnum)
  b=1/meanJsize
  # update lambda which has influence on a,mu,beta
  #diffZV=diff(Z)-diff(V)
  lmbd_pre=100;lmbd_cur=lambda
  if (is.na(lmbd_cur-lmbd_pre))  print("error outside of while loop")
  n=1
  while (abs(lmbd_cur-lmbd_pre)>1e-02) {
    # update m & beta
    VZ=VZupdate(0,V,Z,JJJ,lmbd_cur,delta)
    Q=(diff(VZ[[2]])-diff(VZ[[1]]))/lmbd_cur
    mean1Q=mean(1/Q);meanYQ=mean(Y/Q)
    meanY=mean(Y);meanQ=mean(Q)
    mu=(meanY+meanQ*meanYQ)/(1+meanQ*mean1Q)
    beta=mu*mean1Q-meanYQ
    # update lambda
    deriv_lmbd=lmbdGrad(Y,JJJ,VZ[[1]],Q,mu,beta,lmbd_cur)
    lmbd_pre=lmbd_cur
    lmbd_cur=lmbd_cur+deriv_lmbd/n
    if (lmbd_cur<=0) lmbd_cur=0.01
    #if (lmbd_cur>10) lmbd_cur=10
    n=n+1
    if (is.na(lmbd_cur-lmbd_pre))  print("error in while loop")
  }
  lambda=lmbd_cur
  # estimate a
  a=mean(Jnum)/lambda
  theta=c(lambda,mu,beta,a,b)
  return(theta)
}

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
