#Warwick model

#K = 0.05 1/day, V0 = 0.256 cm^3, rd = 0.2 cm

library(deSolve)

odeFuncW <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    r <- ((3/(4*pi))*x)^(1/3)
    if(r <= rd){
      dx <- K*x
    }
    else{
      dx <- K*x*(1-(1-rd/r)^3)
    }
    return(list(c(dx)))
  })
}

paramW=c(0.05,0.2)
odeparsW = c(K = paramW[1], rd = paramW[2])
odeinitW <- c(x=0.256)


dt=1/24 #hourly
T=10
odeOutW = ode(odeinitW, seq(0,T,length=(T/dt+1)), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
plot(ts(odeOutW[,2],start=0,deltat=dt)) #10 days

#Initial test - simulated data with noise
set.seed(1)
data=odeOutW[,2]+rnorm(T/dt+1,0,0.01) #additive N(0,0.01^2) noise
plot(ts(odeOutW[,2],start=0,deltat=dt)) #10 days
lines(ts(data,start=0,deltat=dt),col=2)
timestest = seq(0,T,length=(T/dt+1))


#Infer K, rd and noise Std Dev

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

loglikeW=function(param=c(0.05,0.2,0.256,0.1),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(K = param[1], rd = param[2])
 odeinit = c(x = param[3])
 out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(y,x,param[4],log=TRUE)))
}

loglikeW(param=c(0.05,0.2,0.256,0.01), y=data, times=timestest)

mhW=function(iters=1000,tune=diag(c(0.01,0.01,0.01,0.01)),init=c(0.3,0.05,0.05,0.4),y=mouse1,times=times1,rt=1e-2,at=1e-3,a=0,b=10)
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
  llikecan=loglikeW(exp(can),y,times,rt,at)
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

matW=mhW(iters=10000, y=data, times = timestest, init=c(0.05,0.2,0.256,0.01))
plot(ts(exp(matW))) #Terrible trace plots from pilot, but results in an ODE solution that fits the data

#update sigtune - use 2000 as burn in
round((1/4)*(2.38^2)*var(matW[2001:10000,]),7)
#           [,1]       [,2]       [,3]       [,4]
#[1,]  0.0236336 -0.0373727  0.0001383  0.0008901
#[2,] -0.0373727  0.0601192 -0.0003594 -0.0012537
#[3,]  0.0001383 -0.0003594  0.0000246 -0.0000440
#[4,]  0.0008901 -0.0012537 -0.0000440  0.0029500


sigtuneW=round((1/4)*(2.38^2)*var(matW[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW[2001:10000,1])) #0.06188379
mean(exp(matW[2001:10000,2])) #0.1405906
mean(exp(matW[2001:10000,3])) #0.2566732
mean(exp(matW[2001:10000,4])) #0.009698501

#re-run
matW=mhW(iters=100000,tune=sigtuneW,init=c(0.0619,0.1406,0.2567,0.0097), y=data,times=timestest) # Acceptance rate 0.116
plot(ts(exp(matW))) #Issues with traces of first two parameters, possible identifiability issues
par(mfrow=c(1,4))
plot(density(exp(matW[,1])))
plot(density(exp(matW[,2])))
plot(density(exp(matW[,3])))
plot(density(exp(matW[,4])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW=apply(exp(matW),2,mean)

odeparsW = c(K = paramMeanW[1], rd = paramMeanW[2])
odeinitW <- c(x = paramMeanW[3])

odeOutW = ode(odeinitW, timestest, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(timestest,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW[4])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW[4])))) #28 days
lines(timestest,odeOutW[,2]+qnorm(0.025,0,paramMeanW[4]),col=2) 
lines(timestest,odeOutW[,2]+qnorm(0.975,0,paramMeanW[4]),col=2) 
lines(timestest,data,type="p") #doesn't fit to data, worse than pilot - gradient too steep

