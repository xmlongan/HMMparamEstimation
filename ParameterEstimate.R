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
      Vlmbd[i+1]=Vlmbd[i+1]+drop((exp(epk*lambda)*epk)%*%Jsize[[i]])
    }
  }
  Qlmbd=-(diff(Vlmbd)+Q)/lambda
  devi=Y-(mu*delta+beta*Q)
  logGlmbd=((-1/2+beta*devi)/Q+(1/2)*(devi/Q)^2)%*%Qlmbd
  return(drop(logGlmbd))
}

ParamUpd <- function(Y,JJJ,V,Z,lambda){
  Jnum=JJJ[[1]];Jsize=JJJ[[2]];Jepoch=JJJ[[3]]
  N=length(Y);delta=1
  # estimate b
  meanJsize=sum(sapply(Jsize,sum))/sum(Jnum)
  b=1/meanJsize
  # update lambda which has influence on a,mu,beta
  #diffZV=diff(Z)-diff(V)
  lmbd.pre=100;lmbd.cur=lambda
  n=1
  while (abs(lmbd.cur-lmbd.pre)>1e-03) {
    # update m & beta
    VZ=VZupdate(0,V,Z,JJJ,lmbd.cur,delta)
    Q=(diff(VZ[[2]])-diff(VZ[[1]]))/lmbd.cur
    mean1Q=mean(1/Q);meanYQ=mean(Y/Q)
    meanY=mean(Y);meanQ=mean(Q)
    mu=(meanY+meanQ*meanYQ)/(1+meanQ*mean1Q)
    beta=mu*mean1Q-meanYQ
    # update lambda
    deriv.lmbd=lmbdGrad(Y,JJJ,VZ[[1]],Q,mu,beta,lmbd.cur)
    lmbd.pre=lmbd.cur
    lmbd.cur=lmbd.cur+deriv.lmbd/n
    # what a bug! since we will use exp(-lambda), so lambda could not be
    # too large!!! Meanwhile lambda should be positive!
    if (lmbd.cur<=0) lmbd.cur=0.01
    if (lmbd.cur>=5) lmbd.cur=5
    n=n+1
  }
  lambda=lmbd.cur
  # estimate a
  a=mean(Jnum)/lambda
  theta=c(lambda,mu,beta,a,b)
  return(theta)
}
