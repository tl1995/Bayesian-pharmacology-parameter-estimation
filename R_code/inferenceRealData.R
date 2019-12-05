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
odeinitE <- c(x=0.001)

times1=c(0,3,7,10,14,17,21,24,28) #in days, starting at day 0: 09/09/2014
odeOutE = ode(odeinitE, times1, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
plot(times1,odeOutE[,2], type="l",xlab="time") #28 days, spaced irregularly

#Data - mouse 1
mouse1=c(0.001,0.001,0.014,0.014,0.115,0.157,0.442,0.561,0.93)


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

loglike=function(param=c(0.05,0.001,0.1),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(K = param[1])
 odeinit = c(x = param[2])
 out2 = ode(odeinit, times, odeFuncE, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(y,x,param[3],log=TRUE)))
}

loglike(param=c(0.05,0.001,1))

mh=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.05,0.001,0.01),y=mouse1,times=times1,rt=1e-2,at=1e-3,a=0,b=10)
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
  llikecan=loglike(exp(can),y,times,rt,at)
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

mat=mh(iters=10000)
plot(ts(exp(mat)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0117940 -0.0416984 -0.0094021
#[2,] -0.0416984  0.1523257  0.0304460
#[3,] -0.0094021  0.0304460  0.1359480
sigtune=round((1/3)*(2.38^2)*var(mat[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat[2001:10000,1])) #0.1397392
mean(exp(mat[2001:10000,2])) #0.0187566
mean(exp(mat[2001:10000,3])) #0.05113891

#re-run
mat=mh(iters=100000,tune=sigtune,init=c(0.1397,0.0188,0.0511))
plot(ts(exp(mat)))
par(mfrow=c(1,3))
plot(density(exp(mat[,1])))
plot(density(exp(mat[,2])))
plot(density(exp(mat[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean=apply(exp(mat),2,mean)

odeparsE = c(K = paramMean[1])
odeinitE <- c(x = paramMean[2])

odeOutE = ode(odeinitE, times1, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times1,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean[3])))) #28 days
lines(times1,odeOutE[,2]+qnorm(0.025,0,paramMean[3]),col=2) 
lines(times1,odeOutE[,2]+qnorm(0.975,0,paramMean[3]),col=2) 
lines(times1,mouse1,type="p")

############################################################################################
## NB - possible evidence of multimodality (or at least very differing convergence times) ##
############################################################################################

##### Mouse 2 #####

mouse2 = c(0.001,0.014,0.034,0.065,0.179,0.335,0.785,1.112,1.792)
times2 = times1 #28 days, spaced irregularly

mat2 = mh(iters=10000, y=mouse2, times=times2)
plot(ts(exp(mat2)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat2[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0094209 -0.0347842  0.0067344
#[2,] -0.0347842  0.1319426 -0.0279234
#[3,]  0.0067344 -0.0279234  0.1886662
sigtune2=round((1/3)*(2.38^2)*var(mat2[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat2[2001:10000,1])) #0.1410833
mean(exp(mat2[2001:10000,2])) #0.03489293
mean(exp(mat2[2001:10000,3])) #0.08250052

#re-run
mat2=mh(iters=100000,tune=sigtune2,init=c(0.1411,0.0349,0.0825), y=mouse2, times=times2)
plot(ts(exp(mat2)))
par(mfrow=c(1,3))
plot(density(exp(mat2[,1])))
plot(density(exp(mat2[,2])))
plot(density(exp(mat2[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean2=apply(exp(mat2),2,mean)

odeparsE = c(K = paramMean2[1])
odeinitE <- c(x = paramMean2[2])

odeOutE = ode(odeinitE, times2, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times2,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean2[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean2[3])))) #28 days
lines(times2,odeOutE[,2]+qnorm(0.025,0,paramMean2[3]),col=2) 
lines(times2,odeOutE[,2]+qnorm(0.975,0,paramMean2[3]),col=2) 
lines(times2,mouse2,type="p")

##### Mouse 3 #####

mouse3 = c(0.001,0.014,0.034,0.037,0.292,0.396,0.802,1.692,3.915)
times3 = times1 #28 days, spaced irregularly

mat3 = mh(iters=10000, y=mouse3, times=times3)
plot(ts(exp(mat3)))

#update sigtune - use 2000 as burn in (NB could maybe use a larger burn in here)
round((1/3)*(2.38^2)*var(mat3[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0146443 -0.0980683  0.0674109
#[2,] -0.0980683  0.6603204 -0.4468973
#[3,]  0.0674109 -0.4468973  0.4347848

sigtune3=round((1/3)*(2.38^2)*var(mat3[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat3[2001:10000,1])) #0.2216144
mean(exp(mat3[2001:10000,2])) #0.008354211
mean(exp(mat3[2001:10000,3])) #0.08638376

#re-run
mat3=mh(iters=100000,tune=sigtune3,init=c(0.2216,0.0084,0.0864), y=mouse3, times=times3)
plot(ts(exp(mat3)))
par(mfrow=c(1,3))
plot(density(exp(mat3[,1])))
plot(density(exp(mat3[,2])))
plot(density(exp(mat3[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean3=apply(exp(mat3),2,mean)

odeparsE = c(K = paramMean3[1])
odeinitE <- c(x = paramMean3[2])

odeOutE = ode(odeinitE, times3, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times3,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean3[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean3[3])))) #28 days
lines(times3,odeOutE[,2]+qnorm(0.025,0,paramMean3[3]),col=2) 
lines(times3,odeOutE[,2]+qnorm(0.975,0,paramMean3[3]),col=2) 
lines(times3,mouse3,type="p")

##### Mouse 4 #####

mouse4 = c(0.014,0.014,0.038,0.049,0.127,0.218,0.841,1.313,2.522)
times4 = times1 #28 days, spaced irregularly

mat4 = mh(iters=10000, y=mouse4, times=times4)
plot(ts(exp(mat4)))

#update sigtune - use 4000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat4[4001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0142765 -0.0649427 -0.0277814
#[2,] -0.0649427  0.2998874  0.1217126
#[3,] -0.0277814  0.1217126  0.2133793

sigtune4=round((1/3)*(2.38^2)*var(mat4[4001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat4[2001:10000,1])) #0.1396548
mean(exp(mat4[2001:10000,2])) #0.1091224
mean(exp(mat4[2001:10000,3])) #0.2630933

#re-run
mat4=mh(iters=100000,tune=sigtune4,init=c(0.1397,0.1091,0.2631), y=mouse4, times=times4)
plot(ts(exp(mat4)))
par(mfrow=c(1,3))
plot(density(exp(mat4[,1])))
plot(density(exp(mat4[,2])))
plot(density(exp(mat4[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean4=apply(exp(mat4),2,mean)

odeparsE = c(K = paramMean4[1])
odeinitE <- c(x = paramMean4[2])

odeOutE = ode(odeinitE, times4, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times4,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean4[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean4[3])))) #28 days
lines(times4,odeOutE[,2]+qnorm(0.025,0,paramMean4[3]),col=2) 
lines(times4,odeOutE[,2]+qnorm(0.975,0,paramMean4[3]),col=2) 
lines(times4,mouse4,type="p")

##### Mouse 5 #####

mouse5 = c(0.014,0.014,0.034,0.014,0.055,0.094,0.144,0.666,1.197)
times5 = times1 #28 days, spaced irregularly

mat5 = mh(iters=10000, y=mouse5, times=times5)
plot(ts(exp(mat5)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat5[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0237366 -0.1536049  0.0106429
#[2,] -0.1536049  1.0061411 -0.0765950
#[3,]  0.0106429 -0.0765950  0.1401290

sigtune5=round((1/3)*(2.38^2)*var(mat5[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat5[2001:10000,1])) #0.2296897
mean(exp(mat5[2001:10000,2])) #0.002260715
mean(exp(mat5[2001:10000,3])) #0.08222394

#re-run
mat5=mh(iters=100000,tune=sigtune5,init=c(0.2297,0.0023,0.0822), y=mouse5, times=times5)
plot(ts(exp(mat5)))
par(mfrow=c(1,3))
plot(density(exp(mat5[,1])))
plot(density(exp(mat5[,2])))
plot(density(exp(mat5[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean5=apply(exp(mat5),2,mean)

odeparsE = c(K = paramMean5[1])
odeinitE <- c(x = paramMean5[2])

odeOutE = ode(odeinitE, times5, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times5,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean5[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean5[3])))) #28 days
lines(times5,odeOutE[,2]+qnorm(0.025,0,paramMean5[3]),col=2) 
lines(times5,odeOutE[,2]+qnorm(0.975,0,paramMean5[3]),col=2) 
lines(times5,mouse5,type="p") #2 points outside of predictive interval

##### Mouse 6 #####

mouse6 = c(0.014,0.13,0.314,0.827,1.987)
times6 = c(0,3,7,10,14) #14 days, spaced irregularly

mat6 = mh(iters=10000, y=mouse6, times=times6)
plot(ts(exp(mat6)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat6[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0237366 -0.1536049  0.0106429
#[2,] -0.1536049  1.0061411 -0.0765950
#[3,]  0.0106429 -0.0765950  0.1401290

sigtune6=round((1/3)*(2.38^2)*var(mat6[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat6[2001:10000,1])) #0.2397417
mean(exp(mat6[2001:10000,2])) #0.06985904
mean(exp(mat6[2001:10000,3])) #0.07945782

#re-run
mat6=mh(iters=100000,tune=sigtune6,init=c(0.2397,0.0699,0.0795), y=mouse6, times=times6)
plot(ts(exp(mat6)))
par(mfrow=c(1,3))
plot(density(exp(mat6[,1])))
plot(density(exp(mat6[,2])))
plot(density(exp(mat6[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean6=apply(exp(mat6),2,mean)

odeparsE = c(K = paramMean6[1])
odeinitE <- c(x = paramMean6[2])

odeOutE = ode(odeinitE, times6, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times6,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean6[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean6[3])))) #14 days
lines(times6,odeOutE[,2]+qnorm(0.025,0,paramMean6[3]),col=2) 
lines(times6,odeOutE[,2]+qnorm(0.975,0,paramMean6[3]),col=2) 
lines(times6,mouse6,type="p")

##### Mouse 7 #####

mouse7 = c(0.014,0.001,0.014,0.001,0.001)
times7 = c(0,3,7,10,14) #14 days, spaced irregularly

mat7 = mh(iters=10000, y=mouse7, times=times7)
plot(ts(exp(mat7)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat7[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  7.3098372  4.3047824 -0.1064047
#[2,]  4.3047824  9.6854595 -0.3074414
#[3,] -0.1064047 -0.3074414  0.2967839


sigtune7=round((1/3)*(2.38^2)*var(mat7[2001:10000,]),7) 

#Initialise at mean values from pilot
mean(exp(mat7[2001:10000,1])) #0.05247587
mean(exp(mat7[2001:10000,2])) #0.001671267
mean(exp(mat7[2001:10000,3])) #0.010041

#re-run
mat7=mh(iters=100000,tune=0.15*sigtune7,init=c(0.0525,0.0011,0.0102), y=mouse7, times=times7) #factor of 0.15 so as to not break the ODE solver
#####Will sometimes still break ODE solver#####
plot(ts(exp(mat7)))
par(mfrow=c(1,3))
plot(density(exp(mat7[,1])))
plot(density(exp(mat7[,2])))
plot(density(exp(mat7[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean7=apply(exp(mat7),2,mean)

odeparsE = c(K = paramMean7[1])
odeinitE <- c(x = paramMean7[2])

odeOutE = ode(odeinitE, times7, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times7,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean7[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean7[3])))) #14 days
lines(times7,odeOutE[,2]+qnorm(0.025,0,paramMean7[3]),col=2) 
lines(times7,odeOutE[,2]+qnorm(0.975,0,paramMean7[3]),col=2) 
lines(times7,mouse7,type="p") #Model poor fit in this case

##### Mouse 8 #####

mouse8 = c(0.014,0.094,0.129,0.369,1.116)
times8 = c(0,3,7,10,14) #14 days, spaced irregularly

mat8 = mh(iters=10000, y=mouse8, times=times8)
plot(ts(exp(mat8)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat8[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0077605 -0.0289077  0.0037414
#[2,] -0.0289077  0.1099964 -0.0164899
#[3,]  0.0037414 -0.0164899  0.1652066


sigtune8=round((1/3)*(2.38^2)*var(mat8[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat8[2001:10000,1])) #0.2736811
mean(exp(mat8[2001:10000,2])) #0.02350815
mean(exp(mat8[2001:10000,3])) #0.03584243

#re-run
mat8=mh(iters=100000,tune=sigtune8,init=c(0.2737,0.0235,0.0358), y=mouse8, times=times8)
plot(ts(exp(mat8)))
par(mfrow=c(1,3))
plot(density(exp(mat8[,1])))
plot(density(exp(mat8[,2])))
plot(density(exp(mat8[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean8=apply(exp(mat8),2,mean)

odeparsE = c(K = paramMean8[1])
odeinitE <- c(x = paramMean8[2])

odeOutE = ode(odeinitE, times8, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times8,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean8[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean8[3])))) #14 days
lines(times8,odeOutE[,2]+qnorm(0.025,0,paramMean8[3]),col=2) 
lines(times8,odeOutE[,2]+qnorm(0.975,0,paramMean8[3]),col=2) 
lines(times8,mouse8,type="p")

##### Mouse 9 #####

mouse9 = c(0.014,0.063,0.2,0.385,0.953)
times9 = c(0,3,7,10,14) #14 days, spaced irregularly

mat9 = mh(iters=10000, y=mouse9, times=times9)
plot(ts(exp(mat9)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat9[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0031314 -0.0099883  0.0031484
#[2,] -0.0099883  0.0325919 -0.0105327
#[3,]  0.0031484 -0.0105327  0.1038667

sigtune9=round((1/3)*(2.38^2)*var(mat9[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat9[2001:10000,1])) #0.2304892
mean(exp(mat9[2001:10000,2])) #0.03655376
mean(exp(mat9[2001:10000,3])) #0.0188087

#re-run
mat9=mh(iters=100000,tune=sigtune9,init=c(0.2305,0.0366,0.0188), y=mouse9, times=times9)
plot(ts(exp(mat9)))
par(mfrow=c(1,3))
plot(density(exp(mat9[,1])))
plot(density(exp(mat9[,2])))
plot(density(exp(mat9[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean9=apply(exp(mat9),2,mean)

odeparsE = c(K = paramMean9[1])
odeinitE <- c(x = paramMean9[2])

odeOutE = ode(odeinitE, times9, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times9,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean9[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean9[3])))) #14 days
lines(times9,odeOutE[,2]+qnorm(0.025,0,paramMean9[3]),col=2) 
lines(times9,odeOutE[,2]+qnorm(0.975,0,paramMean9[3]),col=2) 
lines(times9,mouse9,type="p")

##### Mouse 10 #####

mouse10 = c(0.014,0.069,0.251,0.505,1.104)
times10 = c(0,3,7,10,14) #14 days, spaced irregularly

mat10 = mh(iters=10000, y=mouse10, times=times10)
plot(ts(exp(mat10)))

#update sigtune - use 7000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat10[7001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0073068 -0.0210094 -0.0008161
#[2,] -0.0210094  0.0636717  0.0037411
#[3,] -0.0008161  0.0037411  0.0919575

sigtune10=round((1/3)*(2.38^2)*var(mat10[7001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat10[7001:10000,1])) #0.2136945
mean(exp(mat10[7001:10000,2])) #0.05501696
mean(exp(mat10[7001:10000,3])) #0.04449217

#re-run
mat10=mh(iters=100000,tune=sigtune10,init=c(0.2137,0.0550,0.0445), y=mouse10, times=times10)
plot(ts(exp(mat10)))
par(mfrow=c(1,3))
plot(density(exp(mat10[,1])))
plot(density(exp(mat10[,2])))
plot(density(exp(mat10[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean10=apply(exp(mat10),2,mean)

odeparsE = c(K = paramMean10[1])
odeinitE <- c(x = paramMean10[2])

odeOutE = ode(odeinitE, times10, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times10,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean10[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean10[3])))) #14 days
lines(times10,odeOutE[,2]+qnorm(0.025,0,paramMean10[3]),col=2) 
lines(times10,odeOutE[,2]+qnorm(0.975,0,paramMean10[3]),col=2) 
lines(times10,mouse10,type="p")

##### Mouse 11 #####

mouse11 = c(0.014,0.014,0.014,0.001,0.014,0.001,0.001,0.014,0.034)
times11 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat11 = mh(iters=10000, y=mouse11, times=times11)
plot(ts(exp(mat11)))

#update sigtune - use 8000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat11[8001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0172016 -0.0546784  0.0282341
#[2,] -0.0546784  0.1825194 -0.1198430
#[3,]  0.0282341 -0.1198430  0.5148319

sigtune11=round((1/3)*(2.38^2)*var(mat11[7001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat11[8001:10000,1])) #0.01213732
mean(exp(mat11[8001:10000,2])) #0.009623172
mean(exp(mat11[8001:10000,3])) #0.01165071

#re-run
mat11=mh(iters=100000,tune=0.3*sigtune11,init=c(0.0121,0.0096,0.0117), y=mouse11, times=times11)
#factor of 0.3 here to stop ODE solver breaking - may still break
plot(ts(exp(mat11)))
par(mfrow=c(1,3))
plot(density(exp(mat11[,1])))
plot(density(exp(mat11[,2])))
plot(density(exp(mat11[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean11=apply(exp(mat11),2,mean)

odeparsE = c(K = paramMean11[1])
odeinitE <- c(x = paramMean11[2])

odeOutE = ode(odeinitE, times11, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times11,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean11[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean11[3])))) #28 days
lines(times11,odeOutE[,2]+qnorm(0.025,0,paramMean11[3]),col=2) 
lines(times11,odeOutE[,2]+qnorm(0.975,0,paramMean11[3]),col=2) 
lines(times11,mouse11,type="p") #model doesn't appear to be a great fit

##### Mouse 12 #####

mouse12 = c(0.001,0.014,0.014,0.001,0.001,0.001,0.001,0.014,0.001)
times12 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat12 = mh(iters=10000, y=mouse12, times=times12)
plot(ts(exp(mat12)))

#update sigtune - use 8000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat12[8001:10000,]),7)
#          [,1]       [,2]       [,3]
#[1,] 2.3539390  0.0759080  0.0581788
#[2,] 0.0759080  0.2745429 -0.0436329
#[3,] 0.0581788 -0.0436329  0.1142703

sigtune12=round((1/3)*(2.38^2)*var(mat12[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat12[8001:10000,1])) #0.001922034
mean(exp(mat12[8001:10000,2])) #0.004662528
mean(exp(mat12[8001:10000,3])) #0.006848606

#re-run
mat12=mh(iters=100000,tune=0.2*sigtune12,init=c(0.0019,0.0047,0.0068), y=mouse12, times=times12)
#factor of 0.2 here to stop ODE solver breaking - may still break
plot(ts(exp(mat12)))
par(mfrow=c(1,3))
plot(density(exp(mat12[,1])))
plot(density(exp(mat12[,2])))
plot(density(exp(mat12[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean12=apply(exp(mat12),2,mean)

odeparsE = c(K = paramMean12[1])
odeinitE <- c(x = paramMean12[2])

odeOutE = ode(odeinitE, times12, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times12,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean12[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean12[3])))) #28 days
lines(times12,odeOutE[,2]+qnorm(0.025,0,paramMean12[3]),col=2) 
lines(times12,odeOutE[,2]+qnorm(0.975,0,paramMean12[3]),col=2) 
lines(times12,mouse12,type="p") #model doesn't appear to be a great fit

##### Mouse 13 #####

mouse13 = c(0.001,0.001,0.014,0.001,0.001,0.001,0.001,0.001,0.001)
times13 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat13 = mh(iters=10000, y=mouse13, times=times13)
plot(ts(exp(mat13)))

#update sigtune - use 8000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat13[8001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  1.2613405 -0.4787691 -0.0076811
#[2,] -0.4787691  1.4267334 -0.0139915
#[3,] -0.0076811 -0.0139915  0.1277631

sigtune13=round((1/3)*(2.38^2)*var(mat13[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat13[8001:10000,1])) #0.00000000009180533 - this may cause problems so will initialise at a larger value
mean(exp(mat13[8001:10000,2])) #0.0005001571
mean(exp(mat13[8001:10000,3])) #0.005175234

#re-run
mat13=mh(iters=100000,tune=0.2*sigtune13,init=c(0.0001,0.0005,0.0052), y=mouse13, times=times13)
#factor of 0.2 here to stop ODE solver breaking - may still break
plot(ts(exp(mat13)))
par(mfrow=c(1,3))
plot(density(exp(mat13[,1])))
plot(density(exp(mat13[,2])))
plot(density(exp(mat13[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean13=apply(exp(mat13),2,mean)

odeparsE = c(K = paramMean13[1])
odeinitE <- c(x = paramMean13[2])

odeOutE = ode(odeinitE, times13, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times13,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean13[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean13[3])))) #28 days
lines(times13,odeOutE[,2]+qnorm(0.025,0,paramMean13[3]),col=2) 
lines(times13,odeOutE[,2]+qnorm(0.975,0,paramMean13[3]),col=2) 
lines(times13,mouse13,type="p") #1 point outside of predictive interval, model doesn't appear to be a great fit

##### Mouse 14 #####

mouse14 = c(0.001,0.014,0.014,0.001,0.001,0.001,0.001,0.001,0.014)
times14 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat14 = mh(iters=10000, y=mouse14, times=times14)
plot(ts(exp(mat14)))

#update sigtune - use 8000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat14[8001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.3642371 -0.0336523 -0.0236488
#[2,] -0.0336523  3.2945582  0.0251937
#[3,] -0.0236488  0.0251937  0.1015454

sigtune14=round((1/3)*(2.38^2)*var(mat13[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat14[8001:10000,1])) #0.184321
mean(exp(mat14[8001:10000,2])) #0.00001250664 - this may cause problems with ODE solver so will initialise at a higher value
mean(exp(mat14[8001:10000,3])) #0.00812993

#re-run
mat14=mh(iters=100000,tune=0.1*sigtune12,init=c(0.1843,0.0001,0.0081), y=mouse14, times=times14)
#factor of 0.33 and use of sigtune12 here to stop ODE solver breaking - may still break (NB mouse 12 and mouse 14 have very similar data)
plot(ts(exp(mat14)))
par(mfrow=c(1,3))
plot(density(exp(mat14[,1])))
plot(density(exp(mat14[,2])))
plot(density(exp(mat14[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean14=apply(exp(mat14),2,mean)

odeparsE = c(K = paramMean14[1])
odeinitE <- c(x = paramMean14[2])

odeOutE = ode(odeinitE, times14, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times14,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean14[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean14[3])))) #28 days
lines(times14,odeOutE[,2]+qnorm(0.025,0,paramMean14[3]),col=2) 
lines(times14,odeOutE[,2]+qnorm(0.975,0,paramMean14[3]),col=2) 
lines(times14,mouse14,type="p") #model doesn't appear to be a great fit


##### Mouse 15 #####

mouse15 = c(0.014,0.014,0.014,0.001,0.001,0.001,0.001,0.001,0.001)
times15 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat15 = mh(iters=10000, y=mouse15, times=times15)
plot(ts(exp(mat15)))

#update sigtune - use 8000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat15[8001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  1.7456926  0.4391140 -0.0116827
#[2,]  0.4391140  1.1559962 -0.0131558
#[3,] -0.0116827 -0.0131558  0.0931457

sigtune15=round((1/3)*(2.38^2)*var(mat15[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat15[8001:10000,1])) #0.0001951723
mean(exp(mat15[8001:10000,2])) #0.0000000231629 - this may cause problems with ODE solver so will initialise at a higher value
mean(exp(mat15[8001:10000,3])) #0.008849644

#re-run
mat15=mh(iters=100000,tune=0.33*sigtune15,init=c(0.0002,0.0001,0.0088), y=mouse15, times=times15)
#factor of 0.33 here to stop ODE solver breaking - may still break
plot(ts(exp(mat15)))
par(mfrow=c(1,3))
plot(density(exp(mat15[,1])))
plot(density(exp(mat15[,2])))
plot(density(exp(mat15[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean15=apply(exp(mat15),2,mean)

odeparsE = c(K = paramMean15[1])
odeinitE <- c(x = paramMean15[2])

odeOutE = ode(odeinitE, times15, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times15,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean15[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean15[3])))) #28 days
lines(times15,odeOutE[,2]+qnorm(0.025,0,paramMean15[3]),col=2) 
lines(times15,odeOutE[,2]+qnorm(0.975,0,paramMean15[3]),col=2) 
lines(times15,mouse15,type="p") #model doesn't appear to be a great fit

##### Mouse 16 #####

mouse16 = c(0.014,0.014,0.034,0.014,0.001,0.001,0.001,0.001,0.001)
times16 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat16 = mh(iters=10000, y=mouse16, times=times16)
plot(ts(exp(mat16)))

#update sigtune - use 6000 as burn in (Not sure if it will affect things using a longer burn in?)
round((1/3)*(2.38^2)*var(mat16[6001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  2.5271356 -0.5510065  0.0865627
#[2,] -0.5510065  1.6813768 -0.0872994
#[3,]  0.0865627 -0.0872994  0.1300073

sigtune16=round((1/3)*(2.38^2)*var(mat16[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat15[8001:10000,1])) #0.00207185
mean(exp(mat15[8001:10000,2])) #0.003128534
mean(exp(mat15[8001:10000,3])) #0.00798876

#re-run
mat16=mh(iters=100000,tune=0.25*sigtune16,init=c(0.0021,0.0031,0.0080), y=mouse16, times=times16)
#factor of 0.25 here to stop ODE solver breaking - may still break
plot(ts(exp(mat16)))
par(mfrow=c(1,3))
plot(density(exp(mat16[,1])))
plot(density(exp(mat16[,2])))
plot(density(exp(mat16[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean16=apply(exp(mat16),2,mean)

odeparsE = c(K = paramMean16[1])
odeinitE <- c(x = paramMean16[2])

odeOutE = ode(odeinitE, times16, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times16,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean16[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean16[3])))) #28 days
lines(times16,odeOutE[,2]+qnorm(0.025,0,paramMean16[3]),col=2) 
lines(times16,odeOutE[,2]+qnorm(0.975,0,paramMean16[3]),col=2) 
lines(times16,mouse16,type="p") #1 point outside predictive interval - model doesn't appear to be a great fit


##### Mouse 17 #####

mouse17 = c(0.014,0.014,0.034,0.171,0.328,0.384,0.971,1.513,1.906)
times17 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat17 = mh(iters=10000, y=mouse17, times=times17)
plot(ts(exp(mat17)))

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(mat17[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0394535 -0.1230107 -0.0021871
#[2,] -0.1230107  0.4038682 -0.0130077
#[3,] -0.0021871 -0.0130077  0.1633544

sigtune17=round((1/3)*(2.38^2)*var(mat17[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat17[2001:10000,1])) #0.122259
mean(exp(mat17[2001:10000,2])) #0.07086938
mean(exp(mat17[2001:10000,3])) #0.1674304

#re-run
mat17=mh(iters=100000,tune=sigtune17,init=c(0.1223,0.0709,0.1674), y=mouse17, times=times17)
plot(ts(exp(mat17)))
par(mfrow=c(1,3))
plot(density(exp(mat17[,1])))
plot(density(exp(mat17[,2])))
plot(density(exp(mat17[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean17=apply(exp(mat17),2,mean)

odeparsE = c(K = paramMean17[1])
odeinitE <- c(x = paramMean17[2])

odeOutE = ode(odeinitE, times17, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times17,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean17[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean17[3])))) #28 days
lines(times17,odeOutE[,2]+qnorm(0.025,0,paramMean17[3]),col=2) 
lines(times17,odeOutE[,2]+qnorm(0.975,0,paramMean17[3]),col=2) 
lines(times17,mouse17,type="p")

##### Mouse 18 #####

mouse18 = c(0.014,0.014,0.034,0.039,0.215,0.474,1.05,1.38,1.9)
times18 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat18 = mh(iters=10000, y=mouse18, times=times18)
plot(ts(exp(mat18)))

#update sigtune - use 4000 as burn in
round((1/3)*(2.38^2)*var(mat18[4001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0321457 -0.0929378 -0.0180892
#[2,] -0.0929378  0.2846543  0.0417448
#[3,] -0.0180892  0.0417448  0.1328965

sigtune18=round((1/3)*(2.38^2)*var(mat18[4001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat18[4001:10000,1])) #0.1231565
mean(exp(mat18[4001:10000,2])) #0.06647937
mean(exp(mat18[4001:10000,3])) #0.1716089

#re-run
mat18=mh(iters=100000,tune=sigtune18,init=c(0.1232,0.0665,0.1716), y=mouse18, times=times18)
plot(ts(exp(mat18)))
par(mfrow=c(1,3))
plot(density(exp(mat18[,1])))
plot(density(exp(mat18[,2])))
plot(density(exp(mat18[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean18=apply(exp(mat18),2,mean)

odeparsE = c(K = paramMean18[1])
odeinitE <- c(x = paramMean18[2])

odeOutE = ode(odeinitE, times18, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times18,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean18[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean18[3])))) #28 days
lines(times18,odeOutE[,2]+qnorm(0.025,0,paramMean18[3]),col=2) 
lines(times18,odeOutE[,2]+qnorm(0.975,0,paramMean18[3]),col=2) 
lines(times18,mouse18,type="p")

##### Mouse 19 #####

mouse19 = c(0.055,0.014,0.058,0.034,0.001,0.001,0.001,0.001,0.001)
times19 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat19 = mh(iters=10000, y=mouse19, times=times19)
plot(ts(exp(mat19)))

#update sigtune - use 8000 as burn in
round((1/3)*(2.38^2)*var(mat19[8001:10000,]),7)
#           [,1]      [,2]       [,3]
#[1,]  2.7096302 0.0840409 -0.0015741
#[2,]  0.0840409 0.3650331  0.0336997
#[3,] -0.0015741 0.0336997  0.1184304

sigtune19=round((1/3)*(2.38^2)*var(mat19[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat19[4001:10000,1])) #0.0005943411
mean(exp(mat19[4001:10000,2])) #0.001717956
mean(exp(mat19[4001:10000,3])) #0.02648119

#re-run
mat19=mh(iters=100000,tune=0.3*sigtune19,init=c(0.0006,0.0017,0.0265), y=mouse19, times=times19)
plot(ts(exp(mat19)))
par(mfrow=c(1,3))
plot(density(exp(mat19[,1])))
plot(density(exp(mat19[,2])))
plot(density(exp(mat19[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean19=apply(exp(mat19),2,mean)

odeparsE = c(K = paramMean19[1])
odeinitE <- c(x = paramMean19[2])

odeOutE = ode(odeinitE, times19, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times19,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean19[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean19[3])))) #28 days
lines(times19,odeOutE[,2]+qnorm(0.025,0,paramMean19[3]),col=2) 
lines(times19,odeOutE[,2]+qnorm(0.975,0,paramMean19[3]),col=2) 
lines(times19,mouse19,type="p") #model doesn't appear to be a great fit

##### Mouse 19 #####

mouse19 = c(0.055,0.014,0.058,0.034,0.001,0.001,0.001,0.001,0.001)
times19 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat19 = mh(iters=10000, y=mouse19, times=times19)
plot(ts(exp(mat19)))

#update sigtune - use 8000 as burn in
round((1/3)*(2.38^2)*var(mat19[8001:10000,]),7)
#           [,1]      [,2]       [,3]
#[1,]  2.7096302 0.0840409 -0.0015741
#[2,]  0.0840409 0.3650331  0.0336997
#[3,] -0.0015741 0.0336997  0.1184304

sigtune19=round((1/3)*(2.38^2)*var(mat19[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat19[4001:10000,1])) #0.0005943411
mean(exp(mat19[4001:10000,2])) #0.001717956
mean(exp(mat19[4001:10000,3])) #0.02648119

#re-run
mat19=mh(iters=100000,tune=0.3*sigtune19,init=c(0.0006,0.0017,0.0265), y=mouse19, times=times19)
plot(ts(exp(mat19)))
par(mfrow=c(1,3))
plot(density(exp(mat19[,1])))
plot(density(exp(mat19[,2])))
plot(density(exp(mat19[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean19=apply(exp(mat19),2,mean)

odeparsE = c(K = paramMean19[1])
odeinitE <- c(x = paramMean19[2])

odeOutE = ode(odeinitE, times19, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times19,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean19[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean19[3])))) #28 days
lines(times19,odeOutE[,2]+qnorm(0.025,0,paramMean19[3]),col=2) 
lines(times19,odeOutE[,2]+qnorm(0.975,0,paramMean19[3]),col=2) 
lines(times19,mouse19,type="p") #model doesn't appear to be a great fit

##### Mouse 20 #####

mouse20 = c(0.001,0.014,0.014,0.014,0.001,0.014,0.001,0.001,0.001)
times20 = c(0,3,7,10,14,17,21,24,28) #28 days, spaced irregularly

mat20 = mh(iters=10000, y=mouse20, times=times20)
plot(ts(exp(mat20)))

#update sigtune - use 8000 as burn in
round((1/3)*(2.38^2)*var(mat20[8001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  1.8834356 -0.1517184  0.0203331
#[2,] -0.1517184  0.4054372 -0.0018879
#[3,]  0.0203331 -0.0018879  0.0933780

sigtune20=round((1/3)*(2.38^2)*var(mat20[8001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(mat20[8001:10000,1])) #0.0007496605
mean(exp(mat20[8001:10000,2])) #0.005795502
mean(exp(mat20[8001:10000,3])) #0.006985667

#re-run
mat20=mh(iters=100000,tune=0.3*sigtune20,init=c(0.0007,0.0058,0.007), y=mouse20, times=times20)
plot(ts(exp(mat20)))
par(mfrow=c(1,3))
plot(density(exp(mat20[,1])))
plot(density(exp(mat20[,2])))
plot(density(exp(mat20[,3])))

#predictives?
#Cheat by fixing parameter values at their posterior mean
paramMean20=apply(exp(mat20),2,mean)

odeparsE = c(K = paramMean20[1])
odeinitE <- c(x = paramMean20[2])

odeOutE = ode(odeinitE, times20, odeFuncE, odeparsE, rtol=1e-3, atol=1e-3)
par(mfrow=c(1,1))
plot(times20,odeOutE[,2],xlab="time (days)",type='l',ylim=c(min(odeOutE[,2]+qnorm(0.025,0,paramMean20[3])),max(odeOutE[,2]+qnorm(0.975,0,paramMean20[3])))) #28 days
lines(times20,odeOutE[,2]+qnorm(0.025,0,paramMean20[3]),col=2) 
lines(times20,odeOutE[,2]+qnorm(0.975,0,paramMean20[3]),col=2) 
lines(times20,mouse20,type="p") #model doesn't appear to be a great fit


