#Warwick model

#K = 0.05 1/day, V0 = 0.256 cm^3, rd = 0.2 cm

library(deSolve)

odeFuncW <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    r <- ((3/(4*pi))*x)^(1/3)
    K <- alpha / rd
    if(r <= rd){
      dx <- K*x
    }
    else{
      dx <- K*x*(1-(1-rd/r)^3)
    }
    return(list(c(dx)))
  })
}

paramW=c(0.01, 0.256)
odeinitW <- c(x=paramW[2])
odeparsW = c(alpha = paramW[1], rd = ((3/(4*pi))*paramW[2])^(1/3))



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


#Infer V0, alpha and noise Std Dev

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

loglikeW=function(param=c(0.05,0.256,0.1),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(alpha = param[1], rd = ((3/(4*pi))*param[2])^(1/3))
 odeinit = c(x = param[2])
 out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(y,x,param[3],log=TRUE)))
}

loglikeWMult=function(param=c(0.05,0.256,0.1),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
  odepars = c(alpha = param[1], rd = ((3/(4*pi))*param[2])^(1/3))
  odeinit = c(x = param[2])
  out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
  x=out2[,2] #latent process 
  return(sum(dnorm(log(y),log(x),param[3],log=TRUE)))
}

loglikeW(param=c(0.01,0.256,0.01), y=data, times=timestest)
loglikeWMult(param=c(0.01,0.256,1), y=data, times=timestest)

mhW=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.3,0.05,0.4),y=mouse1,times=times1,rt=1e-2,at=1e-3,a=0,b=10, mult=F)
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
  if(mult==T){
    llikecan=loglikeWMult(exp(can),y,times,rt,at)
  }
  else{
    llikecan=loglikeW(exp(can),y,times,rt,at)
  }
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

matW=mhW(iters=10000, y=data, times = timestest, init=c(0.01,0.256,0.01))
plot(ts(exp(matW))) #Terrible trace plots from pilot, but results in an ODE solution that fits the data

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matW[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0010781 -0.0001541 -0.0004266
#[2,] -0.0001541  0.0000308  0.0000844
#[3,] -0.0004266  0.0000844  0.0033271


sigtuneW=round((1/3)*(2.38^2)*var(matW[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW[2001:10000,1])) #0.009704874
mean(exp(matW[2001:10000,2])) #0.2573074
mean(exp(matW[2001:10000,3])) #0.009637868

#re-run
matW=mhW(iters=100000,tune=sigtuneW,init=c(0.0097,0.2573,0.0096), y=data,times=timestest) # Acceptance rate 0.302
plot(ts(exp(matW))) #Trace plots fine, identifiability issues seem to have gone
par(mfrow=c(1,3))
plot(density(exp(matW[,1])))
plot(density(exp(matW[,2])))
plot(density(exp(matW[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW=apply(exp(matW),2,mean)

odeparsW = c(alpha = paramMeanW[1], rd = ((3/(4*pi))*paramMeanW[2])^(1/3))
odeinitW <- c(x = paramMeanW[2])

odeOutW = ode(odeinitW, timestest, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(timestest,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW[3])))) #28 days
lines(timestest,odeOutW[,2]+qnorm(0.025,0,paramMeanW[3]),col=2) 
lines(timestest,odeOutW[,2]+qnorm(0.975,0,paramMeanW[3]),col=2) 
lines(timestest,data,type="p")

#Data - mouse 1
mouse1=c(0.001,0.001,0.014,0.014,0.115,0.157,0.442,0.561,0.93)
times1=c(0,3,7,10,14,17,21,24,28) #in days, starting at day 0: 09/09/2014

paramW=c(0.01, 0.001)
odeinitW <- c(x=paramW[2])
odeparsW = c(K = paramW[1], rd = ((3/(4*pi))*paramW[2])^(1/3))

matW1=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01))
plot(ts(exp(matW1)))

matW1Mult=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), mult=T)
plot(ts(exp(matW1Mult)))


#predictives from pilot run
paramMeanW1=apply(exp(matW1[2000:10000,]),2,mean)

odeparsW = c(alpha = paramMeanW1[1], rd = ((3/(4*pi))*paramMeanW1[2])^(1/3))
odeinitW <- c(x = paramMeanW1[2])

odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
plot(times1,odeOutW[,2], type="l",xlab="time", ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW1[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW1[3])))) #28 days, spaced irregularly
lines(times1,odeOutW[,2]+qnorm(0.025,0,paramMeanW1[3]),col=2) 
lines(times1,odeOutW[,2]+qnorm(0.975,0,paramMeanW1[3]),col=2) 
lines(times1,mouse1,type="p") #good!

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matW1[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0006697  0.0229576 -0.0004217
#[2,]  0.0229576  1.0697477 -0.0018924
#[3,] -0.0004217 -0.0018924  0.1019723

sigtuneW1=round((1/3)*(2.38^2)*var(matW1[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW1[2001:10000,1])) #0.02654579
mean(exp(matW1[2001:10000,2])) #0.002703788
mean(exp(matW1[2001:10000,3])) #0.03367438


#re-run
matW1=mhW(iters=100000,tune=sigtuneW1, init=c(0.0265, 0.0027, 0.0337)) # Acceptance rate 0.038
plot(ts(exp(matW1)))
par(mfrow=c(1,3))
plot(density(exp(matW1[,1])))
plot(density(exp(matW1[,2])))
plot(density(exp(matW1[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW1=apply(exp(matW1),2,mean)

odeparsW = c(alpha = paramMeanW1[1], rd = ((3/(4*pi))*paramMeanW1[2])^(1/3))
odeinitW <- c(x = paramMeanW1[2])

odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times1,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW1[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW1[3])))) #28 days
lines(times1,odeOutW[,2]+qnorm(0.025,0,paramMeanW1[3]),col = 2) 
lines(times1,odeOutW[,2]+qnorm(0.975,0,paramMeanW1[3]),col = 2) 
lines(times1,mouse1,type="p") #much better than previous attempts, but still fits worse than pilot - gradient now too shallow

#plotting individual solution trajectories
plot(times1, mouse1, type = "n", main = "1000 solutions with plausible parameters", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
odeparsW = c(alpha = exp(matW1[i*100,1]), rd = ((3/(4*pi))*exp(matW1[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matW1)[i*100,2])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(times1, mouse1, type="p")


#Data - mouse 2
mouse2 = c(0.001,0.014,0.034,0.065,0.179,0.335,0.785,1.112,1.792)
times2 = times1 #28 days, starting at day 0: 09/09/2014

paramW=c(0.01, 0.001)
odeinitW <- c(x=paramW[2])
odeparsW = c(K = paramW[1], rd = ((3/(4*pi))*paramW[2])^(1/3))

matW2=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse2, times = times2)
plot(ts(exp(matW2)))

#predictives from pilot run
paramMeanW2=apply(exp(matW2[2000:10000,]),2,mean)

odeparsW = c(alpha = paramMeanW2[1], rd = ((3/(4*pi))*paramMeanW2[2])^(1/3))
odeinitW <- c(x = paramMeanW2[2])

odeOutW = ode(odeinitW, times2, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
plot(times1,odeOutW[,2], type="l",xlab="time", ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW2[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW2[3])))) #28 days, spaced irregularly
lines(times1,odeOutW[,2]+qnorm(0.025,0,paramMeanW2[3]),col=2) 
lines(times1,odeOutW[,2]+qnorm(0.975,0,paramMeanW2[3]),col=2) 
lines(times1,mouse1,type="p") #less good

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matW2[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0003001  0.0095969 -0.0007026
#[2,]  0.0095969  0.3633577 -0.0317889
#[3,] -0.0007026 -0.0317889  0.1354329

sigtuneW2=round((1/3)*(2.38^2)*var(matW2[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW2[2001:10000,1])) #0.03230221
mean(exp(matW2[2001:10000,2])) #0.001624214
mean(exp(matW2[2001:10000,3])) #0.03387839


#re-run
matW2=mhW(iters=100000,tune=sigtuneW2, init=c(0.0323, 0.0016, 0.0339), y = mouse2, times = times2) # Acceptance rate 0.313
plot(ts(exp(matW2)))
par(mfrow=c(1,3))
plot(density(exp(matW2[,1])))
plot(density(exp(matW2[,2])))
plot(density(exp(matW2[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW2=apply(exp(matW2),2,mean)

odeparsW = c(alpha = paramMeanW2[1], rd = ((3/(4*pi))*paramMeanW2[2])^(1/3))
odeinitW <- c(x = paramMeanW2[2])

odeOutW = ode(odeinitW, times2, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times2,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW2[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW2[3])))) #28 days
lines(times2,odeOutW[,2]+qnorm(0.025,0,paramMeanW2[3]),col = 2) 
lines(times2,odeOutW[,2]+qnorm(0.975,0,paramMeanW2[3]),col = 2) 
lines(times2,mouse2,type="p") #much better than previous attempts, but still fits worse than pilot - gradient now too shallow

#plotting individual solution trajectories
plot(times2, mouse2, type = "n")
for(i in 1:1000){
odeparsW = c(alpha = exp(matW2[i*100,1]), rd = ((3/(4*pi))*exp(matW2[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matW2)[i*100,2])
odeOutW = ode(odeinitW, times2, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times2, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}
points(times2, mouse2, type="p")

#Data - mouse 3
mouse3 = c(0.001,0.014,0.034,0.037,0.292,0.396,0.802,1.692,3.915)
times3 = times1 #28 days, starting at day 0: 09/09/2014

matW3=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse3, times = times3)
plot(ts(exp(matW3)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matW3[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0032836  0.0774262 -0.0018747
#[2,]  0.0774262  2.9216238 -0.0534665
#[3,] -0.0018747 -0.0534665  0.1208420


sigtuneW3=round((1/3)*(2.38^2)*var(matW3[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW3[2001:10000,1])) #0.03828567
mean(exp(matW3[2001:10000,2])) #0.001622971
mean(exp(matW3[2001:10000,3])) #0.3376443


#re-run
matW3=mhW(iters=100000,tune=sigtuneW3, init=c(0.0383, 0.0016, 0.3376), y = mouse3, times = times3) # Acceptance rate 0.379
plot(ts(exp(matW3)))
par(mfrow=c(1,3))
plot(density(exp(matW3[,1])))
plot(density(exp(matW3[,2])))
plot(density(exp(matW3[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW3=apply(exp(matW3),2,mean)

odeparsW = c(alpha = paramMeanW3[1], rd = ((3/(4*pi))*paramMeanW3[2])^(1/3))
odeinitW <- c(x = paramMeanW3[2])

odeOutW = ode(odeinitW, times3, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times3,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW3[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW3[3])))) #28 days
lines(times3,odeOutW[,2]+qnorm(0.025,0,paramMeanW3[3]),col = 2) 
lines(times3,odeOutW[,2]+qnorm(0.975,0,paramMeanW3[3]),col = 2) 
lines(times3,mouse3,type="p") #very good fit

#plotting individual solution trajectories
plot(times3, mouse3, type = "n")
for(i in 1:1000){
odeparsW = c(alpha = exp(matW3[i*100,1]), rd = ((3/(4*pi))*exp(matW3[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matW3)[i*100,2])
odeOutW = ode(odeinitW, times3, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times3, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}
points(times3, mouse3, type="p")

#Data - mouse 4
mouse4 = c(0.014,0.014,0.038,0.049,0.127,0.218,0.841,1.313,2.522)
times4 = times1 #28 days, starting at day 0: 09/09/2014

matW4=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse4, times = times4)
plot(ts(exp(matW4)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matW4[2001:10000,]),7)
#          [,1]      [,2]      [,3]
#[1,] 0.0021098 0.0621152 0.0023581
#[2,] 0.0621152 2.3157768 0.0983863
#[3,] 0.0023581 0.0983863 0.1525744


sigtuneW4=round((1/3)*(2.38^2)*var(matW4[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW4[2001:10000,1])) #0.03480289
mean(exp(matW4[2001:10000,2])) #0.002197971
mean(exp(matW4[2001:10000,3])) #0.1420361


#re-run
matW4=mhW(iters=100000,tune=sigtuneW4, init=c(0.0348, 0.0022, 0.1420), y = mouse4, times = times4) # Acceptance rate 0.347
plot(ts(exp(matW4)))
par(mfrow=c(1,3))
plot(density(exp(matW4[,1])))
plot(density(exp(matW4[,2])))
plot(density(exp(matW4[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW4=apply(exp(matW4),2,mean)

odeparsW = c(alpha = paramMeanW4[1], rd = ((3/(4*pi))*paramMeanW4[2])^(1/3))
odeinitW <- c(x = paramMeanW4[2])

odeOutW = ode(odeinitW, times4, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times4,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW4[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW4[3])))) #28 days
lines(times4,odeOutW[,2]+qnorm(0.025,0,paramMeanW4[3]),col = 2) 
lines(times4,odeOutW[,2]+qnorm(0.975,0,paramMeanW4[3]),col = 2) 
lines(times4,mouse4,type="p") #very good fit

#plotting individual solution trajectories
plot(times4, mouse4, type = "n")
for(i in 1:1000){
odeparsW = c(alpha = exp(matW4[i*100,1]), rd = ((3/(4*pi))*exp(matW4[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matW4)[i*100,2])
odeOutW = ode(odeinitW, times4, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times4, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}
points(times4, mouse4, type="p")

#Data - mouse 5
mouse5 = c(0.014,0.014,0.034,0.014,0.055,0.094,0.144,0.666,1.197)
times5 = times1 #28 days, starting at day 0: 09/09/2014

matW5=mhW(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse5, times = times5)
plot(ts(exp(matW5)))

#update sigtune - use 4000 as burn in
round((1/3)*(2.38^2)*var(matW5[4001:10000,]),7)
#           [,1]      [,2]       [,3]
#[1,]  0.0043569 0.0911940 -0.0013990
#[2,]  0.0911940 4.5426636  0.0226334
#[3,] -0.0013990 0.0226334  0.1520813



sigtuneW5=round((1/3)*(2.38^2)*var(matW5[4001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matW5[4001:10000,1])) #0.0242595
mean(exp(matW5[4001:10000,2])) #0.00003696733
mean(exp(matW5[4001:10000,3])) #0.139184


#re-run
matW5=mhW(iters=100000,tune=sigtuneW5, init=c(0.0243, 0.0001, 0.1392), y = mouse5, times = times5) # Acceptance rate 0.347
plot(ts(exp(matW5)))
par(mfrow=c(1,3))
plot(density(exp(matW5[,1])))
plot(density(exp(matW5[,2])))
plot(density(exp(matW5[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMeanW5=apply(exp(matW5),2,mean)

odeparsW = c(alpha = paramMeanW5[1], rd = ((3/(4*pi))*paramMeanW5[2])^(1/3))
odeinitW <- c(x = paramMeanW5[2])

odeOutW = ode(odeinitW, times5, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times5,odeOutW[,2],xlab="time (days)",type='l',ylim=c(min(odeOutW[,2]+qnorm(0.025,0,paramMeanW5[3])),max(odeOutW[,2]+qnorm(0.975,0,paramMeanW5[3])))) #28 days
lines(times5,odeOutW[,2]+qnorm(0.025,0,paramMeanW5[3]),col = 2) 
lines(times5,odeOutW[,2]+qnorm(0.975,0,paramMeanW5[3]),col = 2) 
lines(times5,mouse5,type="p") #very good fit

#plotting individual solution trajectories
plot(times5, mouse5, type = "n", main = "Without random effects")
for(i in 1:1000){
odeparsW = c(alpha = exp(matW5[i*100,1]), rd = ((3/(4*pi))*exp(matW5[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matW5)[i*100,2])
odeOutW = ode(odeinitW, times5, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times5, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}
points(times5, mouse5, type="p")
