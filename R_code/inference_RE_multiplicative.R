#Warwick model

#Without random effects first

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

loglikeWMult=function(param=c(0.05,0.256,0.1),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(alpha = param[1], rd = ((3/(4*pi))*param[2])^(1/3))
 odeinit = c(x = param[2])
 out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(log(y),log(x),param[3],log=TRUE)))
}

loglikeWMult(param=c(0.01,0.256,0.01), y=data, times=timestest)

mhWMult=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.3,0.05,0.4),y=mouse1,times=times1,rt=1e-2,at=1e-3,a=0,b=10)
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
  llikecan=loglikeWMult(exp(can),y,times,rt,at)
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

#Mouse data - defined in inferenceWarwick.R

#mouse 1
matWMult1=mhWMult(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01))
plot(ts(exp(matWMult1)))

#plotting individual solution trajectories
plot(times1, mouse1, type = "n", main = "1000 solutions with plausible parameters", xlab = "Time", ylab = "Volume")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult1[i,1]), rd = ((3/(4*pi))*exp(matWMult1[i,2]))^(1/3))
odeinitW <- c(x = exp(matWMult1)[i,2])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times1,mouse1, type="p")

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matWMult1[2001:10000,]),7)
#          [,1]       [,2]       [,3]
#[1,] 0.0119213  0.0060534  0.0003031
#[2,] 0.0060534  0.4526055 -0.0227306
#[3,] 0.0003031 -0.0227306  0.1602727

sigtuneWMult1 = round((1/3)*(2.38^2)*var(matWMult1[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matWMult1[2001:10000,1])) #0.02440325
mean(exp(matWMult1[2001:10000,2])) #0.0006260203
mean(exp(matWMult1[2001:10000,3])) #0.6141678

#re-run
matWMult1=mhWMult(iters=100000,tune=sigtuneWMult1, init=c(0.0244, 0.0006, 0.6142), y = mouse1, times = times1) # Acceptance rate 0.27552

plot(ts(exp(matWMult1)))
par(mfrow=c(1,3))
plot(density(exp(matWMult1[,1])))
plot(density(exp(matWMult1[,2])))
plot(density(exp(matWMult1[,3])))
par(mfrow=c(1,1))

#plotting individual solution trajectories
plot(times1, mouse1, type = "n", main = "1000 solutions with plausible parameters", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWMult1[i*100,1]), rd = ((3/(4*pi))*exp(matWMult1[i*100,2]))^(1/3))
odeinitW <- c(x = exp(matWMult1)[i*100,2])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times1,mouse1, type="p")

##### Mouse 2 #####

matWMult2=mhWMult(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse2, times = times2) #pilot acc. rate 0.4503
plot(ts(exp(matWMult2)))

#plotting individual solution trajectories
plot(times2, mouse2, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult2[i,1]), rd = ((3/(4*pi))*exp(matWMult2[i,2]))^(1/3))
odeinitW <- c(x = exp(matWMult2)[i,2])
odeOutW = ode(odeinitW, times2, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times2, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times2,mouse2, type="p")

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matWMult2[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0049696 -0.0006777 -0.0012311
#[2,] -0.0006777  0.2681999 -0.0336118
#[3,] -0.0012311 -0.0336118  0.1183607


sigtuneWMult2 = round((1/3)*(2.38^2)*var(matWMult2[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matWMult2[2001:10000,1])) #0.03284564
mean(exp(matWMult2[2001:10000,2])) #0.001566779
mean(exp(matWMult2[2001:10000,3])) #0.3995776

#re-run
matWMult2=mhWMult(iters=100000,tune=sigtuneWMult2, init=c(0.0328, 0.0016, 0.3996), y = mouse2, times = times2) # Acceptance rate 0.30399

plot(ts(exp(matWMult2)))
par(mfrow=c(1,3))
plot(density(exp(matWMult2[,1])))
plot(density(exp(matWMult2[,2])))
plot(density(exp(matWMult2[,3])))
par(mfrow=c(1,1))

#plotting individual solution trajectories
plot(times2, mouse2, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult2[i*10,1]), rd = ((3/(4*pi))*exp(matWMult2[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matWMult2)[i*10,2])
odeOutW = ode(odeinitW, times2, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times2, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times2,mouse2, type="p")

##### Mouse 3 #####

matWMult3=mhWMult(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse3, times = times3) #pilot acc. rate 0.5415
plot(ts(exp(matWMult3)))

#plotting individual solution trajectories
plot(times3, mouse3, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult3[i,1]), rd = ((3/(4*pi))*exp(matWMult3[i,2]))^(1/3))
odeinitW <- c(x = exp(matWMult3)[i,2])
odeOutW = ode(odeinitW, times3, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times3, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times3,mouse3, type="p")

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matWMult3[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0093655 -0.0053472 -0.0007588
#[2,] -0.0053472  0.5358991  0.0065197
#[3,] -0.0007588  0.0065197  0.1180659

sigtuneWMult3 = round((1/3)*(2.38^2)*var(matWMult3[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matWMult3[2001:10000,1])) #0.03505787
mean(exp(matWMult3[2001:10000,2])) #0.001631606
mean(exp(matWMult3[2001:10000,3])) #0.5747707

#re-run
matWMult3=mhWMult(iters=100000,tune=sigtuneWMult3, init=c(0.0351, 0.0016, 0.5748), y = mouse3, times = times3) # Acceptance rate 0.30453

plot(ts(exp(matWMult3)))
par(mfrow=c(1,3))
plot(density(exp(matWMult3[,1])))
plot(density(exp(matWMult3[,2])))
plot(density(exp(matWMult3[,3])))
par(mfrow=c(1,1))

#plotting individual solution trajectories
plot(times3, mouse3, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult3[i*10,1]), rd = ((3/(4*pi))*exp(matWMult3[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matWMult3)[i*10,2])
odeOutW = ode(odeinitW, times3, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times3, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times3,mouse3, type="p")

##### Mouse 4 #####

matWMult4=mhWMult(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse4, times = times4) #pilot acc. rate 0.5613
plot(ts(exp(matWMult4)))

#plotting individual solution trajectories
plot(times4, mouse4, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult4[i,1]), rd = ((3/(4*pi))*exp(matWMult4[i,2]))^(1/3))
odeinitW <- c(x = exp(matWMult4)[i,2])
odeOutW = ode(odeinitW, times4, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times4, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times4,mouse4, type="p")

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matWMult4[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0141958 -0.0096293 -0.0000780
#[2,] -0.0096293  0.2983523  0.0625748
#[3,] -0.0000780  0.0625748  0.1887653

sigtuneWMult4 = round((1/3)*(2.38^2)*var(matWMult4[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matWMult4[2001:10000,1])) #0.033110331
mean(exp(matWMult4[2001:10000,2])) #0.009812673
mean(exp(matWMult4[2001:10000,3])) #0.5394019

#re-run
matWMult4=mhWMult(iters=100000,tune=sigtuneWMult4, init=c(0.0331, 0.098, 0.5394), y = mouse4, times = times4) # Acceptance rate 0.26037

plot(ts(exp(matWMult4)))
par(mfrow=c(1,3))
plot(density(exp(matWMult4[,1])))
plot(density(exp(matWMult4[,2])))
plot(density(exp(matWMult4[,3])))
par(mfrow=c(1,1))

#plotting individual solution trajectories
plot(times4, mouse4, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult4[i*10,1]), rd = ((3/(4*pi))*exp(matWMult4[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matWMult4)[i*10,2])
odeOutW = ode(odeinitW, times4, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times4, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times4,mouse4, type="p")

##### Mouse 5 #####

matWMult5=mhWMult(iters=10000, tune=diag(c(0.01,0.01,0.01)), init = c(0.01,0.001,0.01), y = mouse5, times = times5) #pilot acc. rate 0.5613
plot(ts(exp(matWMult5)))

#plotting individual solution trajectories
plot(times5, mouse5, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult5[i,1]), rd = ((3/(4*pi))*exp(matWMult5[i,2]))^(1/3))
odeinitW <- c(x = exp(matWMult5)[i,2])
odeOutW = ode(odeinitW, times5, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times5, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times5,mouse5, type="p")

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matWMult5[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0417081 -0.0492064 -0.0044257
#[2,] -0.0492064  0.6831208 -0.0046903
#[3,] -0.0044257 -0.0046903  0.1423307

sigtuneWMult5 = round((1/3)*(2.38^2)*var(matWMult5[2001:10000,]),7)

#Initialise at mean values from pilot
mean(exp(matWMult5[2001:10000,1])) #0.02184439
mean(exp(matWMult5[2001:10000,2])) #0.008812296
mean(exp(matWMult5[2001:10000,3])) #0.8149666

#re-run
matWMult5=mhWMult(iters=100000,tune=sigtuneWMult5, init=c(0.0218, 0.0088, 0.815), y = mouse5, times = times5) # Acceptance rate 0.27552

plot(ts(exp(matWMult5)))
par(mfrow=c(1,3))
plot(density(exp(matWMult5[,1])))
plot(density(exp(matWMult5[,2])))
plot(density(exp(matWMult5[,3])))
par(mfrow=c(1,1))

#plotting individual solution trajectories
plot(times5, mouse5, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matWMult5[i*10,1]), rd = ((3/(4*pi))*exp(matWMult5[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matWMult5)[i*10,2])
odeOutW = ode(odeinitW, times5, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times5, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times5,mouse5, type="p")


##### Random effects #####

mhWreMult=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.01,0.256,0.01), hypinit=c(0.02, 0.01),y=data,times=timestest,rt=1e-2,at=1e-3,a=0,b=10, me = 0, d = 0.1, g = 100, h = 1)
{
 ptm = proc.time()
 p=length(init)
 M = dim(y)[2]
 mat=array(0,dim=c(iters,p,M))
 mat[1,,]=log(init)
 alphabar=vector('numeric', length = iters)
 tau=vector('numeric', length = iters)
 alphabar[1] = log(hypinit[1])
 tau[1] = hypinit[2]
 curr=matrix(log(init),nrow=p,ncol=M)
 llikecurr=rep(-1e8, M) #accept first
 count=rep(1, M)
 for(i in 2:iters)
 {
   alpha = mat[(i-1),1,]
   alphabar[i] = rnorm(1, (me*d + M*tau[i-1]*mean(alpha))/(d + M*tau[i-1]), sqrt(1/(d + M*tau[i-1])))
   s2 = 1/M*sum((alpha-mean(alpha))^2)
   tau[i] = rgamma(1, g+M/2, h+M/2*(s2+(mean(alpha)-alphabar[i])^2))
   for(j in 1:M){
    can=rmvn(curr[,j],tune)
    llikecan=loglikeWMult(exp(can),y[,j],times,rt,at)
    laprob=llikecan-llikecurr[j]+logprior(exp(can[1]),alphabar[i],1/sqrt(tau[i]))+logprior(exp(can[2:3]),a,b)-logprior(exp(curr[2:3,j]),a,b)-logprior(exp(curr[1,j]),alphabar[i],1/sqrt(tau[i]))
    if(log(runif(1))<laprob)
    {
     curr[,j]=can
     llikecurr[j]=llikecan
     count[j]=count[j]+1
    }
    mat[i,,j]=curr[,j]
   }
 }
 print(count/iters)
 print(proc.time()-ptm)
 return(list(alphabar = alphabar, tau = tau, mat = mat))
}

matWreMultG1=mhWreMult(iters = 10000,tune=diag(c(0.01,0.01,0.01)),init=c(0.025,0.001,1), hypinit=c(0.025, 0.025),y=mousegroup1,times=times1,rt=1e-2,at=1e-3,a=0,b=10, me = 0, d = 0.1, g = 100, h = 1)

plot(ts(matWreMultG1$alphabar))
plot(density(matWreMultG1$alphabar))
plot(density(exp(matWreMultG1$alphabar)), main = "Multiplicative", xlab=expression(exp(bar(alpha))), ylab="Density")
plot(ts(matWreMultG1$tau))
plot(density(matWreMultG1$tau))
plot(ts(exp(matWreMultG1$mat[,,1])))
plot(ts(exp(matWreMultG1$mat[,,2])))
plot(ts(exp(matWreMultG1$mat[,,3])))
plot(ts(exp(matWreMultG1$mat[,,4])))
plot(ts(exp(matWreMultG1$mat[,,5])))

par(mfrow=c(1,3))
plot(density(exp(matWreMultG1$mat[,1,2])))
plot(density(exp(matWreMultG1$mat[,2,2])))
plot(density(exp(matWreMultG1$mat[,3,2])))
par(mfrow=c(1,1))

round((1/3)*(2.38^2)*var(matWreMultG1$mat[2000:10000,,1]),7)

sigtuneWreMultG1=round((1/3)*(2.38^2)*var(matWreMultG1$mat[2000:10000,,1]),7)

mean(matWreMultG1$alphabar[2000:100000]) #-3.527566
mean(exp(matWreMultG1$alphabar[2000:100000])) #0.02942707
1/var(matWreMultG1$alphabar[2000:100000])  #289.4024
mean(matWreMultG1$tau[2000:100000]) #98.71883
var(matWreMultG1$tau[2000:100000]) #95.15146
mean(exp(matWreMultG1$mat[2000:100000,2,1])) #0.000667888
mean(exp(matWreMultG1$mat[2000:100000,3,1])) #0.5769846

matWreMultG1=mhWreMult(iters = 100000, tune=sigtuneWreMultG1, init=c(0.0294, 0.0005, 0.606), hypinit=c(0.03, 99),y=mousegroup1,times=times1,rt=1e-2,at=1e-3,a=0,b=10, me=-3.48, d=400, g=95, h = 1)

mean(exp(matWreMultG1$mat[6000:100000,1,1]))
mean(exp(matWreMultG1$mat[6000:100000,2,1]))
mean(exp(matWreMultG1$mat[6000:100000,3,1]))

sd(exp(matWreMultG1$mat[6000:100000,1,1]))
sd(exp(matWreMultG1$mat[6000:100000,2,1]))
sd(exp(matWreMultG1$mat[6000:100000,3,1]))

mean(exp(matWreMultG1$mat[6000:100000,1,2]))
mean(exp(matWreMultG1$mat[6000:100000,2,2]))
mean(exp(matWreMultG1$mat[6000:100000,3,2]))

sd(exp(matWreMultG1$mat[6000:100000,1,2]))
sd(exp(matWreMultG1$mat[6000:100000,2,2]))
sd(exp(matWreMultG1$mat[6000:100000,3,2]))

mean(exp(matWreMultG1$mat[6000:100000,1,3]))
mean(exp(matWreMultG1$mat[6000:100000,2,3]))
mean(exp(matWreMultG1$mat[6000:100000,3,3]))

sd(exp(matWreMultG1$mat[6000:100000,1,3]))
sd(exp(matWreMultG1$mat[6000:100000,2,3]))
sd(exp(matWreMultG1$mat[6000:100000,3,3]))

mean(exp(matWreMultG1$mat[6000:100000,1,4]))
mean(exp(matWreMultG1$mat[6000:100000,2,4]))
mean(exp(matWreMultG1$mat[6000:100000,3,4]))

sd(exp(matWreMultG1$mat[6000:100000,1,4]))
sd(exp(matWreMultG1$mat[6000:100000,2,4]))
sd(exp(matWreMultG1$mat[6000:100000,3,4]))

mean(exp(matWreMultG1$mat[6000:100000,1,5]))
mean(exp(matWreMultG1$mat[6000:100000,2,5]))
mean(exp(matWreMultG1$mat[6000:100000,3,5]))

sd(exp(matWreMultG1$mat[6000:100000,1,5]))
sd(exp(matWreMultG1$mat[6000:100000,2,5]))
sd(exp(matWreMultG1$mat[6000:100000,3,5]))

#plotting individual solution trajectories
plot(times1, mouse1, type = "n", main = "Multiplicative", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWreMultG1$mat[i*100,1,1]), rd = ((3/(4*pi))*exp(matWreMultG1$mat[i*100,2,1]))^(1/3))
odeinitW <- c(x = exp(matWreMultG1$mat)[i*100,2,1])
odeOutW = ode(odeinitW, times5, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times5, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(times1, mouse1)

#add mean line
meanodepars2=c(alpha = 0.0263, rd = ((3/(4*pi))*0.000681)^(1/3))
meanodeinit2 <- c(x = 0.000681)
meanodeout2=ode(meanodeinit2,times1,odeFuncW,meanodepars2, rtol=1e-3, atol=1e-3)
lines(times1, meanodeout2[,2], col=2, lty=2)

##### Adding sigma and V0 random effects #####

#NB times must now be a matrix of the same dimension as y

#Also, tune should be an array with the third dimension length equal to the column length of y / times

mhWreMult2=function(iters=1000,tune,init=c(0.01,0.256,0.01), hypinit=rep(0.01,6),y=data,times=timestest,rt=1e-2,at=1e-3, me_a = 0, d_a = 0.1, g_a = 100, h_a = 1, me_s = log(0.5), d_s = 0.1, g_s = 100, h_s = 1, me_V = log(0.001), d_V = 0.1, g_V = 100, h_V = 1)
{
 ptm = proc.time()
 p=length(init)
 M = dim(y)[2]
 mat=array(0,dim=c(iters,p,M))
 mat[1,,]=log(init)
 alphabar=vector('numeric', length = iters)
 tau_alpha=vector('numeric', length = iters)
 V0bar=vector('numeric', length = iters)
 tau_V0=vector('numeric', length = iters)
 sigmabar=vector('numeric', length = iters)
 tau_sig=vector('numeric', length = iters)
 
 alphabar[1] = log(hypinit[1])
 tau_alpha[1] = hypinit[2]
 V0bar[1] = log(hypinit[3])
 tau_V0[1] = hypinit[4]
 sigmabar[1] = log(hypinit[5])
 tau_sig[1] = hypinit[6]
 curr=matrix(log(init),nrow=p,ncol=M)
 llikecurr=rep(-1e8, M) #accept first
 count=rep(1, M)
 for(i in 2:iters)
 {
   alpha = mat[(i-1),1,]
   alphabar[i] = rnorm(1, (me_a*d_a + M*tau_alpha[i-1]*mean(alpha))/(d_a + M*tau_alpha[i-1]), sqrt(1/(d_a + M*tau_alpha[i-1])))
   s2_alpha = 1/M*sum((alpha-mean(alpha))^2)
   tau_alpha[i] = rgamma(1, g_a+M/2, h_a+M/2*(s2_alpha+(mean(alpha)-alphabar[i])^2))
  
   V0 = mat[(i-1),2,]
   V0bar[i] = rnorm(1, (me_V*d_V + M*tau_V0[i-1]*mean(V0))/(d_V + M*tau_V0[i-1]), sqrt(1/(d_V + M*tau_V0[i-1])))
   s2_V0 = 1/M*sum((V0-mean(V0))^2)
   tau_V0[i] = rgamma(1, g_V+M/2, h_V+M/2*(s2_V0+(mean(V0)-V0bar[i])^2))

   sigma = mat[(i-1),3,]
   sigmabar[i] = rnorm(1, (me_s*d_s + M*tau_sig[i-1]*mean(sigma))/(d_s + M*tau_sig[i-1]), sqrt(1/(d_s + M*tau_sig[i-1])))
   s2_sig = 1/M*sum((sigma-mean(sigma))^2)
   tau_sig[i] = rgamma(1, g_s+M/2, h_s+M/2*(s2_sig+(mean(sigma)-sigmabar[i])^2))
 
   for(j in 1:M){
    can=rmvn(curr[,j],tune[,,j])
    llikecan=loglikeWMult(exp(can),na.omit(y[,j]),na.omit(times[,j]),rt,at)
    laprob=llikecan-llikecurr[j]+logprior(exp(can[1]),alphabar[i],1/sqrt(tau_alpha[i]))+logprior(exp(can[2]),V0bar[i],1/sqrt(tau_V0[i]))+logprior(exp(can[3]),sigmabar[i],1/sqrt(tau_sig[i]))-logprior(exp(curr[2,j]),V0bar[i],1/sqrt(tau_V0[i]))-logprior(exp(curr[1,j]),alphabar[i],1/sqrt(tau_alpha[i]))-logprior(exp(curr[3,j]),sigmabar[i],1/sqrt(tau_sig[i]))
    if(log(runif(1))<laprob)
    {
     curr[,j]=can
     llikecurr[j]=llikecan
     count[j]=count[j]+1
    }
    mat[i,,j]=curr[,j]
   }
 }
 print(count/iters)
 print(proc.time()-ptm)
 return(list(alphabar = alphabar, tau_alpha = tau_alpha, V0bar = V0bar, tau_V0 = tau_V0, sigmabar = sigmabar, tau_sigma = tau_sig, mat = mat))
}

mousegroup1=cbind(mouse1,mouse2,mouse3,mouse4,mouse5)
times=cbind(times1,times2,times3,times4,times5)

matWreMultG1_v2=mhWreMult2(iters = 10000,tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,5)),init=c(0.025,0.001,1), hypinit=rep(0.01,6),y=mousegroup1,times=times,rt=1e-2,at=1e-3, me_a = 0, d_a = 0.1, g_a = 100, h_a = 1, me_V = log(0.001), d_V = 0.1, g_V = 100, h_V = 1, me_s = log(0.5), d_s = 0.1, g_s = 100, h_s = 1)

plot(ts(matWreMultG1_v2$alphabar))
plot(density(matWreMultG1_v2$alphabar))
plot(ts(matWreMultG1_v2$tau_alpha))
plot(density(matWreMultG1_v2$tau_alpha))

plot(ts(matWreMultG1_v2$V0bar))
plot(density(matWreMultG1_v2$V0bar))
plot(ts(matWreMultG1_v2$tau_V0))
plot(density(matWreMultG1_v2$tau_V0))

plot(ts(matWreMultG1_v2$sigmabar))
plot(density(matWreMultG1_v2$sigmabar))
plot(ts(matWreMultG1_v2$tau_sigma))
plot(density(matWreMultG1_v2$tau_sigma))


plot(ts(exp(matWreMultG1_v2$mat[,,1])))
plot(ts(exp(matWreMultG1_v2$mat[,,2])))
plot(ts(exp(matWreMultG1_v2$mat[,,3])))
plot(ts(exp(matWreMultG1_v2$mat[,,4])))
plot(ts(exp(matWreMultG1_v2$mat[,,5])))

par(mfrow=c(1,3))
plot(density(exp(matWreMultG1_v2$mat[,1,5])))
plot(density(exp(matWreMultG1_v2$mat[,2,5])))
plot(density(exp(matWreMultG1_v2$mat[,3,5])))
par(mfrow=c(1,1))

round((1/3)*(2.38^2)*var(matWreMultG1_v2$mat[2000:10000,,1]),7)

sigtuneWreMultG1_v2=round((1/3)*(2.38^2)*var(matWreMultG1_v2$mat[2000:10000,,1]),7)

mean(matWreMultG1_v2$alphabar[2000:10000]) #-3.52984
mean(exp(matWreMultG1_v2$alphabar[2000:10000])) #0.02936289
1/var(matWreMultG1_v2$alphabar[2000:10000])  #275.1112
mean(matWreMultG1_v2$tau_alpha[2000:10000]) #99.2185
var(matWreMultG1_v2$tau_alpha[2000:10000]) #100.6878

mean(matWreMultG1_v2$V0bar[2000:10000]) #-5.977753
mean(exp(matWreMultG1_v2$V0bar[2000:10000])) #0.002632109
1/var(matWreMultG1_v2$V0bar[2000:10000])  #13.04945
mean(matWreMultG1_v2$tau_V0[2000:10000]) #99.48241
var(matWreMultG1_v2$tau_V0[2000:10000]) #98.32522

mean(matWreMultG1_v2$sigmabar[2000:10000]) #-0.3635894
mean(exp(matWreMultG1_v2$sigmabar[2000:10000])) #0.7005505
1/var(matWreMultG1_v2$sigmabar[2000:10000])  #64.90713
mean(matWreMultG1_v2$tau_sig[2000:10000]) #99.85751
var(matWreMultG1_v2$tau_sig[2000:10000]) #99.55336

matWreMultG1_v2=mhWreMult2(iters = 100000, tune=sigtuneWreMultG1_v2, init=c(0.0294, 0.0026, 0.7006), hypinit=c(0.029,100,0.003,100,0.701,100),y=mousegroup1,times=times1,rt=1e-2,at=1e-3, me_a = 0, d_a = 0.1, g_a = 100, h_a = 1, me_V = log(0.001), d_V = 0.1, g_V = 100, h_V = 1, me_s = log(0.5), d_s = 0.1, g_s = 100, h_s = 1)

plot(times5, mouse5, type = "n", main = "1000 solutions with plausible parameters", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWreMultG1_v2$mat[i*100,1,5]), rd = ((3/(4*pi))*exp(matWreMultG1_v2$mat[i*100,2,5]))^(1/3))
odeinitW <- c(x = exp(matWreMultG1_v2$mat)[i*100,2,5])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.1))
}
points(times5, mouse5)

#Writing data to file for presentation
write.csv(matWreMultG1$alphabar,file="../data/multiplicativeAlphaBar.csv")

write.csv(matWreMultG1$mat[,,1],file="../data/mouse1outMult.csv")
