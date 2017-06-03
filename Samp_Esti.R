# let's call it Sample and Estimate (SE) contemporarily
rm(list=ls())
# windows
source('D:/学位论文/DissertProject/SampleY_GammaOU.R')
source('D:/学位论文/DissertProject/UpdateState.R')
source('D:/学位论文/DissertProject/ParameterEstimate.R')
source('D:/学位论文/DissertProject/SampleEstimate.R')
# ubuntu
# source('/media/longan/DATA/学位论文/DissertProject/SampleY_GammaOU.R')
# source('/media/longan/DATA/学位论文/DissertProject/UpdateState.R')
# source('/media/longan/DATA/学位论文/DissertProject/ParameterEstimate.R')
# source('/media/longan/DATA/学位论文/DissertProject/SampleEstimate.R')

trueTheta = c(0.2,0,0,1,1);delta=1;N=1000
Y=SampleY(trueTheta,delta,N)

lambda=1;mu=0.5;beta=0.5;a=5;b=5
theta=c(lambda,mu,beta,a,b)

iter=50
thetaIter=theta

for(i in 1:iter){
  VZJ=SmpGammaOU(theta[4],theta[5],theta[1],delta,N)
  # VZJ=list(V,Z,J,Jnum,Jsize,Jepoch)
  theta=SmpEstm(Y,VZJ,theta,delta)you
  thetaIter=rbind(thetaIter,theta)
}
rm(i,iter)
#plot(thetaIter)