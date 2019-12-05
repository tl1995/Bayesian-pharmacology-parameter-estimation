#Exponential model

#K = 0.05 1/day, V0 = 0.256 cm^3

library(deSolve)

odeFuncE <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    dx <- K*x
    return(list(c(dx)))
  })
}

paramE=c(0.05)
odeparsE = c(K = paramE[1])
odeinitE <- c(x=0.256)

dt=1/24 #hourly
T=10
odeOutE = ode(odeinitE, seq(0,T,length=(T/dt+1)), odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
plot(ts(odeOutE[,2],start=0,deltat=dt)) #10 days

#Add noise
set.seed(1)
data=odeOutE[,2]+rnorm(T/dt+1,0,0.01) #additive N(0,0.01^2) noise
plot(ts(odeOutE[,2],start=0,deltat=dt)) #2 days
lines(ts(data,start=0,deltat=dt),col=2)

#Infer K and noise Std Dev

logprior=function(param,a,b)
{
 sum(dnorm(log(param),a,b,log=TRUE)) #N(a,b^2) prior for log params
}

rmvn=function(m,V)
{
  p=length(m)
  z=rnorm(p)
  return(m+t(chol(V))%*%z)
}

loglike=function(param=c(0.05,0.256,0.01),y=data,deltat=1/24,T=10,rt=1e-3,at=1e-3)
{
 odepars = c(K = param[1])
 odeinit = c(x = param[2])
 out2 = ode(odeinit, seq(0,T,length=(T/deltat+1)), odeFuncE, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(y,x,param[3],log=TRUE)))
}

loglike(param=c(0.05,0.256,0.01))

mh=function(iters=1000,tune=diag(c(0.0001,0.0001,0.0001)),init=c(0.05,0.256,0.01),y=data,deltat=1/24,T=10,rt=1e-2,at=1e-3,a=0,b=10)
{
 ptm = proc.time()
 p=length(init)
 mat=matrix(0,ncol=p,nrow=iters)
 mat[1,]=log(init)
 curr=log(init)
 llikecurr=-1e8 #accept first
 count=1
 for(i in 2:iters)
 {
  can=rmvn(curr,tune)
  llikecan=loglike(exp(can),y,deltat,T,rt,at)
  laprob=llikecan-llikecurr+logprior(exp(can),a,b)-logprior(exp(curr),a,b)
  if(log(runif(1))<laprob)
  {
   curr=can
   llikecurr=llikecan
   count=count+1
  }
  mat[i,]=curr
 }
 print(proc.time()-ptm)
 print(count/iters) 
 return(mat)
}

mat=mh(iters=10000,tune=diag(c(0.0001,0.0001,0.0001)))
plot(ts(mat))
#update sigtune
round((1/3)*(2.38^2)*var(mat),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0003558 -0.0001018  0.0000439
#[2,] -0.0001018  0.0000357 -0.0000174
#[3,]  0.0000439 -0.0000174  0.0036081
sigtune=round((1/3)*(2.38^2)*var(mat),7)
#re-run
mat=mh(iters=10000,tune=sigtune)
plot(ts(mat))
par(mfrow=c(1,3))
plot(density(exp(mat[,1])))
lines(c(0.05,0.05),c(0,1000))
plot(density(exp(mat[,2])))
lines(c(0.256,0.256),c(0,1000))
plot(density(exp(mat[,3])))
lines(c(0.01,0.01),c(0,1000))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean=apply(exp(mat),2,mean)

odeparsE = c(K = paramMean[1])
odeinitE <- c(x = paramMean[2])

odeOutE = ode(odeinitE, seq(0,T,length=(T/dt+1)), odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
plot(ts(odeOutE[,2],start=0,deltat=dt)) #10 days
lines(ts(odeOutE[,2]+qnorm(0.025,0,paramMean[3]),start=0,deltat=dt)) 
lines(ts(odeOutE[,2]+qnorm(0.975,0,paramMean[3]),start=0,deltat=dt)) 
lines(seq(0,T,length=(T/dt+1)),data,type="p")


