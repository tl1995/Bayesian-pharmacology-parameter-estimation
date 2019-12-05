#Code for reproducing graphs from MCMC output
#from both applications

#Check working directory is the plot_reproduction_code folder
#change if necessary!
getwd()

library(deSolve)
#Helper function used in ODE solver when plotting solutions
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

par(mfrow = c(1,1))
#log(alpha) illustrative prior
plot(seq(-3,3,length.out=100),dnorm(seq(-3,3,length.out=100)),type="l",axes=F,xlab="",ylab="")
axis(1,at=c(-3,-1.8,0.2,1,3), labels=c("",expression(alpha[1],alpha[2],alpha[3]),""))
points(c(-1.8,0.2,1),dnorm(c(-1.8,0.2,1)),pch="x",col=2)
lines(c(-1.8,-1.8,NA,0.2,0.2,NA,1,1),c(0,dnorm(-1.8),NA,0,dnorm(0.2),NA,0,dnorm(1)),lty=3)

dev.copy(pdf,"../TL_presentation/prior_illustration1.pdf")
dev.off()

#alphabar and tau illustrative priors
plot(seq(-5,5,length.out=100),dnorm(seq(-5,5,length.out=100)),type="l",xaxt="n",yaxt="n",main="",xlab="",ylab="")
axis(1,at=0,labels=expression(b),cex.axis=2,padj=1)

dev.copy(pdf,"../TL_presentation/prior_illustration2.pdf")
dev.off()

plot(seq(0,20,length.out=100),dgamma(seq(0,20,length.out=100),shape=2,rate=1/2),type="l",xaxt="n",yaxt="n",main="",xlab="",ylab="")
axis(1,at=0,labels=expression(0),cex.axis=2,padj=1)

dev.copy(pdf,"../TL_presentation/prior_illustration3.pdf")
dev.off()

#Prior densities application 1
#V0 and sigma
plot(x=seq(-40,40),y=dnorm(seq(-40,40),sd=10),main=expression("Priors for" ~ v[0] ~ "and" ~ sigma),type="l",ylim=c(0,0.1),xlab="Parameter",ylab="Density",las=1)

dev.copy(pdf,"../TL_presentation/vagueprior.pdf")
dev.off()

#alphabar
plot(x=seq(-3.5,-3.45,length.out=1000),y=dnorm(seq(-3.5,-3.45,length.out=1000),mean=-3.48,sd=1/400),main=expression("Prior for" ~ bar(alpha)),type="l",xlab=expression(bar(alpha)),ylab="Density",las=1)

dev.copy(pdf,"../TL_presentation/alphabarprior.pdf")
dev.off()

#tau
plot(x=seq(50,150),y=dgamma(seq(50,150),shape=95,rate=1),main=expression("Prior for" ~ tau),type="l",xlab=expression(tau),ylab="Density",ylim=c(0,0.1),las=1)

dev.copy(pdf,"../TL_presentation/tauprior.pdf")
dev.off()


#alphabar posterior density plots
additiveAlphaBar= read.csv("../data/additiveAlphaBar.csv",header=T)
multiplicativeAlphaBar= read.csv("../data/multiplicativeAlphaBar.csv",header=T)
d1 = density(exp(as.vector(additiveAlphaBar$x)))
plot_p(d1$x, d1$y, main = "Posterior distribution", xlab="", ylab="Density", type = "l")
legend("topright", legend=c("Additive"), col=1,lty=1, text.col = 1,
       bty = "n")

dev.copy(postscript,"../TL_presentation/alphaBarDensities1.eps")
dev.off()

plot_p(d1$x, d1$y, main = "Posterior distribution", xlab="", ylab="Density",
       type = "l")
d2 = density(exp(as.vector(multiplicativeAlphaBar$x)))
lines(d2$x, d2$y, col=2, bty = "n")

legend("topright", legend=c("Additive","Multiplicative"), col=c(1,2),lty=1, bty = "n",
       text.col = 1:2 )

dev.copy(postscript,"../TL_presentation/alphaBarDensities2.eps")
dev.off()

#mouse1 solution plots
mouse1outAdd=read.csv("../data/mouse1out.csv", header=T)
mouse1outMult=read.csv("../data/mouse1outMult.csv", header=T)

mouse1=c(0.001,0.001,0.014,0.014,0.115,0.157,0.442,0.561,0.930)
times1=c(0,3,7,10,14,17,21,24,28)

#Graphing application 1 results - alpha
alphameansadd=c(0.0257,0.0295,0.0346,0.0323,0.0269)
alphasdadd=c(0.0013,0.0017,0.0022,0.0019,0.0015)

alphameansmult=c(0.0263,0.0322,0.0334,0.0309,0.0268)
alphasdmult=c(0.0018,0.0016,0.0021,0.0020,0.0025)

plot(alphameansadd,c(15,12,9,6,3),yaxt="n",xlim=c(0.02,0.04),ylim=c(0,15),pch=4,main=expression("Posterior means for" ~ alpha), xlab="Mean", ylab="")
arrows(x0=alphameansadd-alphasdadd,y0=c(15,12,9,6,3),x1=alphameansadd+alphasdadd,y1=c(15,12,9,6,3),code=3,length=0.05,angle=90)
points(alphameansmult,c(14,11,8,5,2),pch=4,col=2)
arrows(x0=alphameansmult-alphasdmult,y0=c(14,11,8,5,2),x1=alphameansmult+alphasdmult,y1=c(14,11,8,5,2),code=3,length=0.05,angle=90,col=2)
axis(2,at=c(14.5,11.5,8.5,5.5,2.5),labels=expression(alpha[1],alpha[2],alpha[3],alpha[4],alpha[5]),las=1)
legend("topright",legend=c("Additive","Multiplicative"),col=c(1,2),pch=4)

dev.copy(pdf,"../TL_presentation/alphaResults.pdf")
dev.off()

#Graphing application 1 results - V0
V0meansadd=c(0.001815,0.000454,0.000595,0.000576,0.002529)
V0sdadd=c(0.0023,0.0012,0.0056,0.0018,0.0051)

V0meansmult=c(0.000681,0.001763,0.001662,0.008850,0.007688)
V0sdmult=c(0.0005,0.0007,0.0013,0.0037,0.0050)

plot(V0meansadd,c(15,12,9,6,3),yaxt="n",xlim=c(-0.01,0.015),ylim=c(0,15),pch=4,main=expression("Posterior means for" ~ V[0]), xlab="Mean", ylab="")
arrows(x0=V0meansadd-V0sdadd,y0=c(15,12,9,6,3),x1=V0meansadd+V0sdadd,y1=c(15,12,9,6,3),code=3,length=0.05,angle=90)
points(V0meansmult,c(14,11,8,5,2),pch=4,col=2)
arrows(x0=V0meansmult-V0sdmult,y0=c(14,11,8,5,2),x1=V0meansmult+V0sdmult,y1=c(14,11,8,5,2),code=3,length=0.05,angle=90,col=2)
axis(2,at=c(14.5,11.5,8.5,5.5,2.5),labels=expression(V[0]^1,V[0]^2,V[0]^3,V[0]^4,V[0]^5),las=1)
legend("topright",legend=c("Additive","Multiplicative"),col=c(1,2),pch=4)

dev.copy(pdf,"../TL_presentation/V0Results.pdf")
dev.off()

#Graphing application 1 results - sigma
sigmeansadd=c(0.0337,0.0310,0.3827,0.1479,0.1363)
sigsdadd=c(0.0102,0.0092,0.1218,0.0445,0.0418)

sigmeansmult=c(0.6311,0.4157,0.5951,0.5371,0.9117)
sigsdmult=c(0.2192,0.1316,0.1971,0.1641,0.2974)

plot(sigmeansadd,c(15,12,9,6,3),yaxt="n",xlim=c(0,1.2),ylim=c(0,15),pch=4,main=expression("Posterior means for" ~ sigma), xlab="Mean", ylab="")
arrows(x0=sigmeansadd-sigsdadd,y0=c(15,12,9,6,3),x1=sigmeansadd+sigsdadd,y1=c(15,12,9,6,3),code=3,length=0.05,angle=90)
points(sigmeansmult,c(14,11,8,5,2),pch=4,col=2)
arrows(x0=sigmeansmult-sigsdmult,y0=c(14,11,8,5,2),x1=sigmeansmult+sigsdmult,y1=c(14,11,8,5,2),code=3,length=0.05,angle=90,col=2)
axis(2,at=c(14.5,11.5,8.5,5.5,2.5),labels=expression(sigma[1],sigma[2],sigma[3],sigma[4],sigma[5]),las=1)
legend("topright",legend=c("Additive","Multiplicative"),col=c(1,2),pch=4)

dev.copy(pdf,"../TL_presentation/sigmaResults.pdf")
dev.off()

#Posterior illustrations mouse 1
#alpha
plot(density(exp(mouse1outAdd$V1)),main=expression("Posterior density for" ~ alpha[1]),xlab=expression(alpha[1]), las=1)
points(c(exp(mouse1outAdd$V1[100])),0,pch="x",col=2, cex=2)

dev.copy(pdf,"../TL_presentation/alphaposts1.pdf")
dev.off()

points(c(exp(mouse1outAdd$V1[200])),0,pch="x",col=4,cex=2)

dev.copy(pdf,"../TL_presentation/alphaposts2.pdf")
dev.off()

points(c(exp(mouse1outAdd$V1[30000])),0,pch="x",col=6,cex=2)

dev.copy(pdf,"../TL_presentation/alphaposts3.pdf")
dev.off()

#V0
plot(density(exp(mouse1outAdd$V2)),main=expression("Posterior density for" ~ V[0]^1),xlab=expression(V[0]^1), las=1)
points(c(exp(mouse1outAdd$V2[100])),0,pch="x",col=2,cex=2)

dev.copy(pdf,"../TL_presentation/V0posts1.pdf")
dev.off()

points(c(exp(mouse1outAdd$V2[200])),0,pch="x",col=4,cex=2)

dev.copy(pdf,"../TL_presentation/V0posts2.pdf")
dev.off()

points(c(exp(mouse1outAdd$V2[30000])),0,pch="x",col=6,cex=2)

dev.copy(pdf,"../TL_presentation/V0posts3.pdf")
dev.off()


#Additive solutions 1
plot(times1, mouse1, type = "n", main = "Solutions", xlab = "Time", ylab = "Volume", las=1)
odeparsW = c(alpha = exp(mouse1outAdd$V1[100]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[100]))^(1/3))
odeinitW <- c(x = exp(mouse1outAdd$V2)[100])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2])
points(times1, mouse1)


dev.copy(pdf,"../TL_presentation/mouse1sols1.pdf")
dev.off()

#Additive solutions 2
odeparsW = c(alpha = exp(mouse1outAdd$V1[200]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[200]))^(1/3))
odeinitW <- c(x = exp(mouse1outAdd$V2)[200])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2])

dev.copy(pdf,"../TL_presentation/mouse1sols2.pdf")
dev.off()

#Additive solutions 3
odeparsW = c(alpha = exp(mouse1outAdd$V1[30000]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[30000]))^(1/3))
odeinitW <- c(x = exp(mouse1outAdd$V2)[30000])
odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times1, odeOutW[,2])

dev.copy(pdf,"../TL_presentation/mouse1sols3.pdf")
dev.off()

par(mfrow=c(1,2))
#Additive + multiplicative solutions 1:10
plot(times1, mouse1, type = "n", main = "Additive", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:10){
  odeparsW = c(alpha = exp(mouse1outAdd$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outAdd$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2],col=gray(0.8,alpha=0.8))
}
points(times1, mouse1)

plot(times1, mouse1, type = "n", main = "Multiplicative", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:10){
  odeparsW = c(alpha = exp(mouse1outMult$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outMult$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outMult$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2], col=gray(0.8,alpha=0.8))
}
points(times1, mouse1)

dev.copy(pdf,"../TL_presentation/mouse1sols10.pdf")
dev.off()

#Additive + multiplicative solutions 1:100
plot(times1, mouse1, type = "n", main = "Additive", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:100){
  odeparsW = c(alpha = exp(mouse1outAdd$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outAdd$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2],col=gray(0.8,alpha=0.5))
}
points(times1, mouse1)

plot(times1, mouse1, type = "n", main = "Multiplicative", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:100){
  odeparsW = c(alpha = exp(mouse1outMult$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outMult$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outMult$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2], col=gray(0.8,alpha=0.5))
}
points(times1, mouse1)

dev.copy(pdf,"../TL_presentation/mouse1sols100.pdf")
dev.off()

#Additive + multiplicative solutions 1:1000
plot(times1, mouse1, type = "n", main = "Additive", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:1000){
  odeparsW = c(alpha = exp(mouse1outAdd$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outAdd$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2],col=gray(0.8,alpha=0.2))
}
points(times1, mouse1)

plot(times1, mouse1, type = "n", main = "Multiplicative", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:1000){
  odeparsW = c(alpha = exp(mouse1outMult$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outMult$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outMult$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2], col=gray(0.8,alpha=0.2))
}
points(times1, mouse1)

dev.copy(pdf,"../TL_presentation/mouse1sols1000.pdf")
dev.off()

#Additive solutions
plot_p(times1, mouse1, type = "n", main = "Additive", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:1000){
  odeparsW = c(alpha = exp(mouse1outAdd$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outAdd$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outAdd$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(times1, mouse1)

#Adding mean line
#means can be found using eg mean(exp(mouse1outAdd$V1))
meanodepars1=c(alpha = 0.0257, rd = ((3/(4*pi))*0.001815)^(1/3))
meanodeinit1 <- c(x = 0.001815)
meanodeout1=ode(meanodeinit1,times1,odeFuncW,meanodepars1, rtol=1e-3, atol=1e-3)
lines(times1, meanodeout1[,2], col=2, lty=2)

#Multiplicative solutions
plot_p(times1, mouse1, type = "n", main = "Multiplicative", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
  odeparsW = c(alpha = exp(mouse1outMult$V1[i*100]), rd = ((3/(4*pi))*exp(mouse1outMult$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(mouse1outMult$V2)[i*100])
  odeOutW = ode(odeinitW, times1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(times1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(times1, mouse1)

#Adding mean line
#means can be found using eg mean(exp(mouse1outMult$V1))
meanodepars2=c(alpha = 0.0263, rd = ((3/(4*pi))*0.000681)^(1/3))
meanodeinit2 <- c(x = 0.000681)
meanodeout2=ode(meanodeinit2,times1,odeFuncW,meanodepars2, rtol=1e-3, atol=1e-3)
lines(times1, meanodeout2[,2], col=2, lty=2)
par(mfrow=c(1,1))

dev.copy(pdf,"../TL_presentation/SolsMouse1withMean.pdf")
dev.off()

#alphabar posterior densities for Gao data
gcAlphaBar=read.csv("../data/gcAlphaBar.csv", header=T)
crcAlphaBar=read.csv("../data/crcAlphaBar.csv", header=T)
brcaAlphaBar=read.csv("../data/brcaAlphaBar.csv", header=T)
pdacAlphaBar=read.csv("../data/pdacAlphaBar.csv", header=T)
nsclcAlphaBar=read.csv("../data/nsclcAlphaBar.csv", header=T)
cmAlphaBar=read.csv("../data/cmAlphaBar.csv", header=T)

plot(density(gcAlphaBar$x), xlim=c(min(gcAlphaBar$x), max(nsclcAlphaBar$x)), ylim = c(0, 5), main = expression("Posterior densities of the " ~ bar(alpha)), xlab=expression(bar(alpha)), las=1)
lines(density(crcAlphaBar$x), col = 2)
lines(density(brcaAlphaBar$x), col = 3)
lines(density(pdacAlphaBar$x), col = 4)
lines(density(nsclcAlphaBar$x), col = 5)
lines(density(cmAlphaBar$x), col = 6)

legend("topleft", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

dev.copy(pdf,"../TL_presentation/alphaBarPosts.pdf")
dev.off()

#V0bar posterior densities for Gao data
gcV0Bar=read.csv("../data/gcV0Bar.csv", header=T)
crcV0Bar=read.csv("../data/crcV0Bar.csv", header=T)
brcaV0Bar=read.csv("../data/brcaV0Bar.csv", header=T)
pdacV0Bar=read.csv("../data/pdacV0Bar.csv", header=T)
nsclcV0Bar=read.csv("../data/nsclcV0Bar.csv", header=T)
cmV0Bar=read.csv("../data/cmV0Bar.csv", header=T)

plot(density(gcV0Bar$x), xlim=c(min(nsclcV0Bar$x), max(gcV0Bar$x)), ylim = c(0, 16), main = expression("Posterior densities of the " ~ bar(V[0])), xlab=expression(bar(V[0])),las=1)
lines(density(crcV0Bar$x), col = 2)
lines(density(brcaV0Bar$x), col = 3)
lines(density(pdacV0Bar$x), col = 4)
lines(density(nsclcV0Bar$x), col = 5)
lines(density(cmV0Bar$x), col = 6)

legend("topleft", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

dev.copy(pdf,"../TL_presentation/V0BarPosts.pdf")
dev.off()

#sigmabar posterior densities for Gao data
gcSigmaBar=read.csv("../data/gcSigmaBar.csv", header=T)
crcSigmaBar=read.csv("../data/crcSigmaBar.csv", header=T)
brcaSigmaBar=read.csv("../data/brcaSigmaBar.csv", header=T)
pdacSigmaBar=read.csv("../data/pdacSigmaBar.csv", header=T)
nsclcSigmaBar=read.csv("../data/nsclcSigmaBar.csv", header=T)
cmSigmaBar=read.csv("../data/cmSigmaBar.csv", header=T)

plot(density(gcSigmaBar$x), xlim=c(min(brcaSigmaBar$x), -1), ylim = c(0, 12), main = expression("Posterior densities of the " ~ bar(sigma)), xlab=expression(bar(sigma)),las=1)
lines(density(crcSigmaBar$x), col = 2)
lines(density(brcaSigmaBar$x), col = 3)
lines(density(pdacSigmaBar$x), col = 4)
lines(density(nsclcSigmaBar$x), col = 5)
lines(density(cmSigmaBar$x), col = 6)

legend("topright", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

dev.copy(pdf,"../TL_presentation/sigmaBarPosts.pdf")
dev.off()

#Alternate "error bar" plots of posteriors
#alphabar
alphabars=cbind(gcAlphaBar$x,crcAlphaBar$x,pdacAlphaBar$x,brcaAlphaBar$x,nsclcAlphaBar$x,cmAlphaBar$x)

plot_p(apply(alphabars,2,mean),c(6,5,4,3,2,1),yaxt="n",xlim=c(-2.2,-1),ylim=c(0,7),pch=4,main=expression("Posterior means for" ~ bar(alpha)), xlab="Mean", ylab="")
arrows(x0=apply(alphabars,2,mean)-apply(alphabars,2,sd),y0=c(6,5,4,3,2,1),x1=apply(alphabars,2,mean)+apply(alphabars,2,sd),y1=c(6,5,4,3,2,1),code=3,length=0.05,angle=90)
axis(2,at=c(6,5,4,3,2,1),labels=c("GC", "CRC", "PDAC", "BRCA", "NSCLC", "CM"),las=1)

dev.copy(pdf, "../TL_presentation/alphaBarMeans.pdf")
dev.off()

#V0
V0bars=cbind(gcV0Bar$x,crcV0Bar$x,pdacV0Bar$x,brcaV0Bar$x,nsclcV0Bar$x,cmV0Bar$x)

plot(apply(V0bars,2,mean),c(6,5,4,3,2,1),yaxt="n",xlim=c(5.2,5.8),ylim=c(0,6),pch=4,main=expression("Posterior means for" ~ bar(V[0])), xlab="Mean", ylab="")
arrows(x0=apply(V0bars,2,mean)-apply(V0bars,2,sd),y0=c(6,5,4,3,2,1),x1=apply(V0bars,2,mean)+apply(V0bars,2,sd),y1=c(6,5,4,3,2,1),code=3,length=0.05,angle=90)
axis(2,at=c(6,5,4,3,2,1),labels=c("GC", "CRC", "PDAC", "BRCA", "NSCLC", "CM"),las=1)

dev.copy(pdf, "../TL_presentation/V0BarMeans.pdf")
dev.off()

#sigma
sigmabars=cbind(gcSigmaBar$x,crcSigmaBar$x,pdacSigmaBar$x,brcaSigmaBar$x,nsclcSigmaBar$x,cmSigmaBar$x)

plot(apply(sigmabars,2,mean),c(6,5,4,3,2,1),yaxt="n",xlim=c(-2.4,-1.6),ylim=c(0,6),pch=4,main=expression("Posterior means for" ~ bar(sigma)), xlab="Mean", ylab="")
arrows(x0=apply(sigmabars,2,mean)-apply(sigmabars,2,sd),y0=c(6,5,4,3,2,1),x1=apply(sigmabars,2,mean)+apply(sigmabars,2,sd),y1=c(6,5,4,3,2,1),code=3,length=0.05,angle=90)
axis(2,at=c(6,5,4,3,2,1),labels=c("GC", "CRC", "PDAC", "BRCA", "NSCLC", "CM"),las=1)

dev.copy(pdf, "../TL_presentation/sigmaBarMeans.pdf")
dev.off()


#Solution plots for first mouse
#requires reading in and extracting the data
library(readxl)

nm = read_xlsx("../data/nm.3954-S2.xlsx", sheet=4)
nmu = nm[nm$Treatment== 'untreated',]
nmu = nmu[,c(2,3,4,6)]

nmuGC = nmu[nmu[,1]=="GC",]

l1GC=which(nmuGC[,4]==0)
l2GC=which(nmuGC[,4]==0)[-1]

ldiffGC=l2GC-l1GC[-length(l1GC)]
max(ldiffGC) #28

dfdataGC = matrix(NA, nrow=28,ncol = 45)
dftimeGC = matrix(NA, nrow=28,ncol = 45)

for(i in 1:44){
  dfdataGC[1:ldiffGC[i],i] = as.matrix(nmuGC[l1GC[i]:(l1GC[(i+1)]-1),3])
  dftimeGC[1:ldiffGC[i],i] = as.matrix(nmuGC[l1GC[i]:(l1GC[(i+1)]-1),4])
}

dfdataGC[1:8,45] = as.matrix(nmuGC[495:502,3])
dftimeGC[1:8,45] = as.matrix(nmuGC[495:502,4])

#Now read in MCMC output
gcOut1=read.csv("../data/gcOut1.csv",header=T)

#Plotting individual solutions
plot_p(dftimeGC[,1], dfdataGC[,1], type = "n", main = "GC tumour", xlab = "Time", ylab = "Volume", las=1)
for(i in 1:1000){
  odeparsW = c(alpha = exp(gcOut1$V1[i*100]), rd = ((3/(4*pi))*exp(gcOut1$V2[i*100]))^(1/3))
  odeinitW <- c(x = exp(gcOut1$V2)[i*100])
  odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(dftimeGC[,1], dfdataGC[,1])

#Bayesian means
mean(exp(gcOut1$V1)) #0.2119637
mean(exp(gcOut1$V2)) #217.865
odeparsW = c(alpha = 0.2119637, rd = ((3/(4*pi))*217.865)^(1/3))
odeinitW <- c(x = 217.865)
odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = 2, lty=2)

#Frequentist means
#Can be found from the first entry in the .xlsm file in data folder
odeparsW = c(alpha = 3.051384362, rd = 52.403884)
odeinitW <- c(x = 209.4769908909)
odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = 1, lty=3)

legend("topleft", legend = c("Bayesian estimate", "Frequentist estimate"), col=c(2,1), lty=c(2,3), cex=0.75)

dev.copy(pdf,"../TL_presentation/GC_sols.pdf")
dev.off()
