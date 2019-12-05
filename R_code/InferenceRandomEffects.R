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
data=matrix(nrow=T/dt+1, ncol=5)
set.seed(1)
for(i in 1:5){
data[,i] = odeOutW[,2]+rnorm(T/dt+1,0,0.005) #additive N(0,0.01^2) noise
}
plot(ts(odeOutW[,2],start=0,deltat=dt)) #10 days
for(i in 1:5){
lines(ts(data[,i],start=0,deltat=dt),col=i)
}
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

loglikeW=function(param=c(0.01,0.256,0.01),y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(alpha = param[1], rd = ((3/(4*pi))*param[2])^(1/3))
 odeinit = c(x = param[2])
 out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(y,x,param[3],log=TRUE)))
}

loglikeW(param=c(0.01,0.256,0.01), y=data[,1], times=timestest)

##Initial algorithm: loop over all mice to perform inference for all

mhWre=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.01,0.256,0.01),y=data,times=timestest,rt=1e-2,at=1e-3,a=0,b=10)
{
 ptm = proc.time()
 p=length(init)
 M = dim(y)[2]
 mat=array(0,dim=c(iters,p,M))
 mat[1,,]=log(init)
 curr=matrix(log(init),nrow=p,ncol=M)
 llikecurr=rep(-1e8, M) #accept first
 count=rep(1, M)
 for(i in 2:iters)
 {
 for(j in 1:M){
  can=rmvn(curr[,j],tune)
  llikecan=loglikeW(exp(can),y[,j],times,rt,at)
  laprob=llikecan-llikecurr[j]+logprior(exp(can),a,b)-logprior(exp(curr[,j]),a,b)
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
 return(mat)
}

#NB mousegroup1 defined later in code
matWre=mhWre(iters=10000, y=mousegroup1, times = times1, init=c(0.01,0.001,0.01))
plot(ts(exp(matWre[,,1])))

#plotting individual solution trajectories
plot(times1, mouse1, type = "p")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWre[i*10,1,1]), rd = ((3/(4*pi))*exp(matWre[i*10,2,1]))^(1/3))
odeinitW <- c(x = exp(matWre)[i*10,2,1])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}

round((1/3)*(2.38^2)*var(matWre[2001:10000,,5]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0039813 -0.0026377 -0.0038540
#[2,] -0.0026377  2.2799877  0.1890040
#[3,] -0.0038540  0.1890040  0.1382287

sigtuneWre=round((1/3)*(2.38^2)*var(matWre[2001:10000,,5]),7)

exp(mean(matWre[2001:10000,1,5])) #0.02740125
exp(mean(matWre[2001:10000,2,5])) #0.003462436
exp(mean(matWre[2001:10000,3,5])) #0.1364861

matWre=mhWre(iters=10000, y=mousegroup1, times = times1, init=c(0.0274,0.0035,0.0136), tune = sigtuneWre)
plot(ts(exp(matWre[,,1])))
par(mfrow=c(1,3))
plot(density(exp(matWre[,1,1])))
plot(density(exp(matWre[,2,1])))
plot(density(exp(matWre[,3,1])))

#plotting individual solution trajectories
plot(times5, mouse5, type = "p")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWre[i*10,1,5]), rd = ((3/(4*pi))*exp(matWre[i*10,2,5]))^(1/3))
odeinitW <- c(x = exp(matWre)[i*10,2,5])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.25))
}

### Gibbs sampling step for FCDs

gibbsWre = function(iters = 1000, alpha = c(0.01,0.015,0.02,0.025,0.03), init = c(0.02,0.01), b = 0.02,d = 0.1,g = 0.01,h = 1){
  p = length(init)
  M = length(alpha)
  mat = matrix(0, nrow=iters, ncol = p)
  thetabar=log(init[1])
  tau=init[2]
  mat[1,] = c(thetabar, tau)
  for(i in 2:iters){
    thetabar = rnorm(1, (b*d + M*tau*mean(log(alpha)))/(d + M*tau), sqrt(1/(d + M*tau)))
    s2 = 1/M*sum((log(alpha)-mean(log(alpha)))^2)
    tau = rgamma(1, g+M/2, h+M/2*(s2+(mean(log(alpha))-thetabar)^2))
    mat[i,] = c(thetabar, tau)
  }
  return(mat)
}

gibbsmat = gibbsWre()
plot(ts(gibbsmat))
plot(density(gibbsmat[,1]))
plot(density(gibbsmat[,2]))

###NEXT STEP: put gibbs sampling and mh algorithm together to hopefully create the full random effects algorithm, test on simulated data

mhWre2=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.01,0.256,0.01), hypinit=c(0.02, 0.01),y=data,times=timestest,rt=1e-2,at=1e-3,a=0,b=10, me = 0, d = 0.1, g = 100, h = 1)
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
    llikecan=loglikeW(exp(can),y[,j],times,rt,at)
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

matWreTest=mhWre2()
plot(ts(matWreTest$alphabar))
plot(ts(matWreTest$tau))
plot(ts(matWreTest$mat[,,1]))
plot(ts(matWreTest$mat[,,2]))
plot(ts(matWreTest$mat[,,3]))
plot(ts(matWreTest$mat[,,4]))
plot(ts(matWreTest$mat[,,5]))

#Real data
mousegroup1=cbind(mouse1,mouse2,mouse3,mouse4,mouse5)

matWreGroup1=mhWre2(iters = 10000,tune=diag(c(0.01,0.01,0.01)),init=c(0.025,0.001,0.04), hypinit=c(0.025, 0.025),y=mousegroup1,times=times1,rt=1e-2,at=1e-3,a=0,b=10, me = -3.48, d = 400, g = 95, h = 1)

plot(ts(matWreGroup1$alphabar))
plot(density(matWreGroup1$alphabar))
plot(density(exp(matWreGroup1$alphabar)), main = "Additive", xlab=expression(exp(bar(alpha))), ylab="Density")
plot(ts(matWreGroup1$tau))
plot(density(matWreGroup1$tau))
plot(ts(exp(matWreGroup1$mat[,,1])))
plot(ts(exp(matWreGroup1$mat[,,2])))
plot(ts(exp(matWreGroup1$mat[,,3])))
plot(ts(exp(matWreGroup1$mat[,,4])))
plot(ts(exp(matWreGroup1$mat[,,5])))

par(mfrow=c(1,3))
plot(density(exp(matWreGroup1$mat[,1,5])))
plot(density(exp(matWreGroup1$mat[,2,5])))
plot(density(exp(matWreGroup1$mat[,3,5])))
par(mfrow=c(1,1))

round((1/3)*(2.38^2)*var(matWreGroup1$mat[6000:10000,,1]),7)
round((1/3)*(2.38^2)*var(matWreGroup1$mat[6000:10000,,2]),7) #very different tuning matrices, based on two different individuals, could present a problem

sigtuneWre1=round((1/3)*(2.38^2)*var(matWreGroup1$mat[6000:10000,,4]),7)

mean(matWreGroup1$alphabar[6000:10000]) #-3.480828
exp(mean(matWreGroup1$alphabar[6000:10000])) #0.03078192
1/var(matWreGroup1$alphabar[6000:10000])  #409.224
mean(matWreGroup1$tau[6000:10000]) #95.34823
var(matWreGroup1$tau[6000:10000]) #90.24714
exp(mean(matWreGroup1$mat[6000:10000,2,4])) #0.001457407
exp(mean(matWreGroup1$mat[6000:10000,3,4])) #0.1249147



matWreGroup1=mhWre2(iters = 100000, tune=sigtuneWre1, init=c(0.0298, 0.00005, 0.1562), hypinit=c(0.031, 95),y=mousegroup1,times=times1,rt=1e-2,at=1e-3,a=0,b=10, me=-3.48, d=400, g=95, h = 1)

mean(exp(matWreGroup1$mat[6000:100000,1,1]))
mean(exp(matWreGroup1$mat[6000:100000,2,1]))
mean(exp(matWreGroup1$mat[6000:100000,3,1]))

sd(exp(matWreGroup1$mat[6000:100000,1,1]))
sd(exp(matWreGroup1$mat[6000:100000,2,1]))
sd(exp(matWreGroup1$mat[6000:100000,3,1]))

mean(exp(matWreGroup1$mat[6000:100000,1,2]))
mean(exp(matWreGroup1$mat[6000:100000,2,2]))
mean(exp(matWreGroup1$mat[6000:100000,3,2]))

sd(exp(matWreGroup1$mat[6000:100000,1,2]))
sd(exp(matWreGroup1$mat[6000:100000,2,2]))
sd(exp(matWreGroup1$mat[6000:100000,3,2]))

mean(exp(matWreGroup1$mat[6000:100000,1,3]))
mean(exp(matWreGroup1$mat[6000:100000,2,3]))
mean(exp(matWreGroup1$mat[6000:100000,3,3]))

sd(exp(matWreGroup1$mat[6000:100000,1,3]))
sd(exp(matWreGroup1$mat[6000:100000,2,3]))
sd(exp(matWreGroup1$mat[6000:100000,3,3]))

mean(exp(matWreGroup1$mat[6000:100000,1,4]))
mean(exp(matWreGroup1$mat[6000:100000,2,4]))
mean(exp(matWreGroup1$mat[6000:100000,3,4]))

sd(exp(matWreGroup1$mat[6000:100000,1,4]))
sd(exp(matWreGroup1$mat[6000:100000,2,4]))
sd(exp(matWreGroup1$mat[6000:100000,3,4]))

mean(exp(matWreGroup1$mat[6000:100000,1,5]))
mean(exp(matWreGroup1$mat[6000:100000,2,5]))
mean(exp(matWreGroup1$mat[6000:100000,3,5]))

sd(exp(matWreGroup1$mat[6000:100000,1,5]))
sd(exp(matWreGroup1$mat[6000:100000,2,5]))
sd(exp(matWreGroup1$mat[6000:100000,3,5]))

#plotting individual solution trajectories
plot(times1, mouse1, type = "n", main = "Additive", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
odeparsW = c(alpha = exp(matWreGroup1$mat[i*100,1,1]), rd = ((3/(4*pi))*exp(matWreGroup1$mat[i*100,2,1]))^(1/3))
odeinitW <- c(x = exp(matWreGroup1$mat)[i*100,2,1])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(times1, mouse1)

#add mean line
meanodepars1=c(alpha = 0.0257, rd = ((3/(4*pi))*0.001815)^(1/3))
meanodeinit1 <- c(x = 0.001815)
meanodeout1=ode(meanodeinit1,times1,odeFuncW,meanodepars1, rtol=1e-3, atol=1e-3)
lines(times1, meanodeout1[,2], col=2, lty=2)

#Writing data to file for presentation
#write.table(matWreGroup1$alphabar,file="../data/additiveAlphaBar.txt")
#advise against using above, below seems better

write.csv(matWreGroup1$alphabar,file="../data/additiveAlphaBar.csv")

write.csv(matWreGroup1$mat[,,1],file="../data/mouse1out.csv")
