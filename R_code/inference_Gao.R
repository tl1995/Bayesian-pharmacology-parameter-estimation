#Requires functions from InferenceRandomEffects.R and inference_RE_multiplicative.R

GC1 = c(250.4,281.0,313.1,344.1,305.0,447.2,559.9,862.1,1115.5,1712.6)
times_GC1 = c(0,3,6,9,12,15,19,24,28,33)

plot(times_GC1, GC1)


#Parameters from the spreadsheet
paramW=c(0.05822821, 52.403884, 209.476991)
odeinitW <- c(x=paramW[3])
odeparsW = c(alpha = paramW[1]*paramW[2], rd = paramW[2])

odeOutW = ode(odeinitW, times_GC1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times_GC1, odeOutW[,2]) #49 days

matGC1=mhW(iters=10000, y=GC1, times = times_GC1, init=c(0.31,250.4,1))
#Pilot acc. rate 0.1971
plot(ts(exp(matGC1))) 

#plotting individual solution trajectories
plot(times_GC1, GC1, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matGC1[i,1]), rd = ((3/(4*pi))*exp(matGC1[i,2]))^(1/3))
odeinitW <- c(x = exp(matGC1)[i,2])
odeOutW = ode(odeinitW, times_GC1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times_GC1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times_GC1, GC1, type = "p") #Actually looks ok

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matGC1[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0036332 -0.0153224  0.0002290
#[2,] -0.0153224  0.0922569 -0.0146302
#[3,]  0.0002290 -0.0146302  0.1234147

sigtuneGC1 = round((1/3)*(2.38^2)*var(matGC1[2001:10000,]),7)

matGC1=mhW(iters=100000, tune = sigtuneGC1, y=GC1, times = times_GC1, init=c(0.31,250.4,1)) #acc. rate 0.29337

plot(ts(exp(matGC1))) #good looking trace plots

mean(exp(matGC1[20001:100000,1])) #0.2506245
sd(exp(matGC1[20001:100000,1])) #0.010426

mean(exp(matGC1[20001:100000,2])) #152.4618
sd(exp(matGC1[20001:100000,2])) #27.8922

mean(exp(matGC1[20001:100000,3])) #86.61721
sd(exp(matGC1[20001:100000,3])) #26.06714

#plotting individual solution trajectories
plot(times_GC1, GC1, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matGC1[i*10,1]), rd = ((3/(4*pi))*exp(matGC1[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matGC1)[i*10,2])
odeOutW = ode(odeinitW, times_GC1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times_GC1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times_GC1, GC1, type = "p") #Actually looks ok

##### Multiplicative Gaussian noise #####

matGC1Mult=mhWMult(iters=10000, y=GC1, times = times_GC1, init=c(0.31,250.4,1))
#Pilot acc. rate 0.409
plot(ts(exp(matGC1Mult))) 

#plotting individual solution trajectories
plot(times_GC1, GC1, type = "n", main = "Without random effects")
for(i in 1:1000){
odeparsW = c(alpha = exp(matGC1Mult[i*10,1]), rd = ((3/(4*pi))*exp(matGC1Mult[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matGC1Mult)[i*10,2])
odeOutW = ode(odeinitW, times_GC1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times_GC1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times_GC1, GC1, type = "p") #Actually looks ok

#update sigtune - use 2000 as burn in
round((1/3)*(2.38^2)*var(matGC1Mult[2001:10000,]),7)
#           [,1]       [,2]       [,3]
#[1,]  0.0130404 -0.0104095 -0.0006631
#[2,] -0.0104095  0.0228162  0.0009064
#[3,] -0.0006631  0.0009064  0.1323522


sigtuneGC1Mult = round((1/3)*(2.38^2)*var(matGC1Mult[2001:10000,]),7)

matGC1Mult=mhWMult(iters=100000, tune = sigtuneGC1Mult, y=GC1, times = times_GC1, init=c(0.31,250.4,1)) #acc. rate 0.30861

plot(ts(exp(matGC1Mult))) #good looking trace plots

mean(exp(matGC1Mult[20001:100000,1])) #0.2157504
sd(exp(matGC1Mult[20001:100000,1])) #0.01824197

mean(exp(matGC1Mult[20001:100000,2])) #211.8286
sd(exp(matGC1Mult[20001:100000,2])) #24.37914

mean(exp(matGC1Mult[20001:100000,3])) #0.1926777
sd(exp(matGC1Mult[20001:100000,3])) #0.05652216

#plotting individual solution trajectories
plot(times_GC1, GC1, type = "n", main = "Without random effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(matGC1Mult[i*10,1]), rd = ((3/(4*pi))*exp(matGC1Mult[i*10,2]))^(1/3))
odeinitW <- c(x = exp(matGC1Mult)[i*10,2])
odeOutW = ode(odeinitW, times_GC1, odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(times_GC1, odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(times_GC1, GC1, type = "p") #Actually looks ok

###### Performing inference on separate types of tumours ######

library(readxl)

nm = read_xlsx("~/epsrc_az/data/nm.3954-S2.xlsx", sheet=4)
nmu = nm[nm$Treatment== 'untreated',]
nmu = nmu[,c(2,3,4,6)]

#different types of tumours
table(nmu[,1])

l1=which(nmu[,4]==0)
l2=which(nmu[,4]==0)[-1]

ldiff=l2-l1[-length(l1)]
max(ldiff) #38

dfdata = matrix(NA, nrow=38,ncol = 229)
dftime = matrix(NA, nrow=38,ncol = 229)

for(i in 1:228){
dfdata[1:ldiff[i],i] = as.matrix(nmu[l1[i]:(l1[(i+1)]-1),3])
dftime[1:ldiff[i],i] = as.matrix(nmu[l1[i]:(l1[(i+1)]-1),4])
}

dfdata[1:5,229] = as.matrix(nmu[2555:2559,3])
dftime[1:5,229] = as.matrix(nmu[2555:2559,4])

test1=mhWMult(iters=10000, y=na.omit(dfdata[,1]), times = na.omit(dftime[,1]), init=c(0.31,250.4,1))

gaoOut = list()

for(i in 1:229){
gaoOut[[i]] = mhWMult(iters=10000, y=na.omit(dfdata[,i]), times = na.omit(dftime[,i]), init=c(0.31,250.4,1))
}

sigtunelist = list()

for(i in 1:229){
sigtunelist[[i]] = round((1/3)*(2.38^2)*var(gaoOut[[i]][2001:10000,]),7)
}

for(i in 1:229){
gaoOut[[i]] = mhWMult(iters=100000, tune = sigtunelist[[i]], y=na.omit(dfdata[,i]), times = na.omit(dftime[,i]), init=c(0.31,250.4,1))
print(i)
}

#Random effects on GC tumors

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

outGC=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,45)), y=dfdataGC, times = dftimeGC, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outGC$alphabar))
plot(density(outGC$alphabar))
mean(exp(outGC$alphabar)) #0.1296671
plot(ts(outGC$tau_alpha))
plot(density(outGC$tau_alpha))
mean(outGC$tau_alpha) #2.326653

plot(ts(outGC$V0bar))
plot(density(outGC$V0bar))
exp(mean(outGC$V0bar)) #259.7564
plot(ts(outGC$tau_V0))
plot(density(outGC$tau_V0))
mean(outGC$tau_V0) #21.4168

plot(ts(outGC$sigmabar))
plot(density(outGC$sigmabar))
exp(mean(outGC$sigmabar)) #0.1685822
plot(ts(outGC$tau_sigma))
plot(density(outGC$tau_sigma))
mean(outGC$tau_sigma) #866.1494

plot(ts(exp(outGC$mat[,,1])))

sigtuneGC = array(0, dim = c(3,3,dim(dfdataGC)[2]))

for(i in 1:dim(dfdataGC)[2]){
sigtuneGC[,,i]=round((1/3)*(2.38^2)*var(outGC$mat[2001:10000,,i]),7)
}

outGC=mhWreMult2(iters=100000, tune = sigtuneGC, y=dfdataGC, times = dftimeGC, init=c(0.31,250.4,1), hypinit = c(0.13, 2.33, 260, 21, 0.17, 870), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Creating matrix of CRC tumor data

nmuCRC = nmu[nmu[,1]=="CRC",]

l1CRC=which(nmuCRC[,4]==0)
l2CRC=which(nmuCRC[,4]==0)[-1]

ldiffCRC=l2CRC-l1CRC[-length(l1CRC)]
max(ldiffCRC) #28

dfdataCRC = matrix(NA, nrow=28,ncol = 45)
dftimeCRC = matrix(NA, nrow=28,ncol = 45)

for(i in 1:44){
dfdataCRC[1:ldiffCRC[i],i] = as.matrix(nmuCRC[l1CRC[i]:(l1CRC[(i+1)]-1),3])
dftimeCRC[1:ldiffCRC[i],i] = as.matrix(nmuCRC[l1CRC[i]:(l1CRC[(i+1)]-1),4])
}

dfdataCRC[1:8,45] = as.matrix(nmuCRC[525:532,3])
dftimeCRC[1:8,45] = as.matrix(nmuCRC[525:532,4])

#random effects on CRC data

outCRC=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,45)), y=dfdataCRC, times = dftimeCRC, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outCRC$alphabar))
plot(density(outCRC$alphabar))
exp(mean(outCRC$alphabar)) #0.1465525
plot(ts(outCRC$tau_alpha))
plot(density(outCRC$tau_alpha))
mean(outCRC$tau_alpha) #2.618819

plot(ts(outCRC$V0bar))
plot(density(outCRC$V0bar))
exp(mean(outCRC$V0bar)) #240.9422
plot(ts(outCRC$tau_V0))
plot(density(outCRC$tau_V0))
mean(outCRC$tau_V0) #44.93012

plot(ts(outCRC$sigmabar))
plot(density(outCRC$sigmabar))
exp(mean(outCRC$sigmabar)) #0.1193981
plot(ts(outCRC$tau_sigma))
plot(density(outCRC$tau_sigma))
mean(outCRC$tau_sigma) #222.9695

plot(ts(exp(outCRC$mat[,,1])))

sigtuneCRC = array(0, dim = c(3,3,dim(dfdataCRC)[2]))

for(i in 1:dim(dfdataCRC)[2]){
sigtuneCRC[,,i]=round((1/3)*(2.38^2)*var(outCRC$mat[2001:10000,,i]),7)
}

outCRC=mhWreMult2(iters=100000, tune = sigtuneCRC, y=dfdataCRC, times = dftimeCRC, init=c(0.31,250.4,1), hypinit = c(0.13, 2.33, 260, 21, 0.17, 870), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Creating matrix of BRCA tumor data

nmuBRCA = nmu[nmu[,1]=="BRCA",]

l1BRCA=which(nmuBRCA[,4]==0)
l2BRCA=which(nmuBRCA[,4]==0)[-1]

length(l1BRCA) #39

ldiffBRCA=l2BRCA-l1BRCA[-length(l1BRCA)]
max(ldiffBRCA) #38

dfdataBRCA = matrix(NA, nrow=38,ncol = 39)
dftimeBRCA = matrix(NA, nrow=38,ncol = 39)

for(i in 1:38){
dfdataBRCA[1:ldiffBRCA[i],i] = as.matrix(nmuBRCA[l1BRCA[i]:(l1BRCA[(i+1)]-1),3])
dftimeBRCA[1:ldiffBRCA[i],i] = as.matrix(nmuBRCA[l1BRCA[i]:(l1BRCA[(i+1)]-1),4])
}

dfdataBRCA[1:5,39] = as.matrix(nmuBRCA[485:489,3])
dftimeBRCA[1:5,39] = as.matrix(nmuBRCA[485:489,4])

#random effects on BRCA data

outBRCA=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,39)), y=dfdataBRCA, times = dftimeBRCA, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outBRCA$alphabar))
plot(density(outBRCA$alphabar))
exp(mean(outBRCA$alphabar)) #0.1743601
plot(ts(outBRCA$tau_alpha))
plot(density(outBRCA$tau_alpha))
mean(outBRCA$tau_alpha) #2.25312

plot(ts(outBRCA$V0bar))
plot(density(outBRCA$V0bar))
exp(mean(outBRCA$V0bar)) #253.0061
plot(ts(outBRCA$tau_V0))
plot(density(outBRCA$tau_V0))
mean(outBRCA$tau_V0) #48.55369

plot(ts(outBRCA$sigmabar))
plot(density(outBRCA$sigmabar))
exp(mean(outBRCA$sigmabar)) #0.1055054
plot(ts(outBRCA$tau_sigma))
plot(density(outBRCA$tau_sigma))
mean(outBRCA$tau_sigma) #869.446

plot(ts(exp(outBRCA$mat[,,1])))

sigtuneBRCA = array(0, dim = c(3,3,dim(dfdataBRCA)[2]))

for(i in 1:dim(dfdataBRCA)[2]){
sigtuneBRCA[,,i]=round((1/3)*(2.38^2)*var(outBRCA$mat[2001:10000,,i]),7)
}

outBRCA=mhWreMult2(iters=100000, tune = sigtuneBRCA, y=dfdataBRCA, times = dftimeBRCA, init=c(0.31,250.4,1), hypinit = c(0.17, 2.25, 250, 48, 0.11, 870), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Creating matrix of CM tumor data

nmuCM = nmu[nmu[,1]=="CM",]

l1CM=which(nmuCM[,4]==0)
l2CM=which(nmuCM[,4]==0)[-1]

length(l1CM) #33
l1CM

ldiffCM=l2CM-l1CM[-length(l1CM)]
max(ldiffCM) #20

dfdataCM = matrix(NA, nrow=20,ncol = 33)
dftimeCM = matrix(NA, nrow=20,ncol = 33)

for(i in 1:32){
dfdataCM[1:ldiffCM[i],i] = as.matrix(nmuCM[l1CM[i]:(l1CM[(i+1)]-1),3])
dftimeCM[1:ldiffCM[i],i] = as.matrix(nmuCM[l1CM[i]:(l1CM[(i+1)]-1),4])
}

dfdataCM[1:6,33] = as.matrix(nmuCM[263:268,3])
dftimeCM[1:6,33] = as.matrix(nmuCM[263:268,4])

#random effects on CM data

outCM=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,33)), y=dfdataCM, times = dftimeCM, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outCM$alphabar))
plot(density(outCM$alphabar))
exp(mean(outCM$alphabar)) #0.2758023
plot(ts(outCM$tau_alpha))
plot(density(outCM$tau_alpha))
mean(outCM$tau_alpha) #5.365993

plot(ts(outCM$V0bar))
plot(density(outCM$V0bar))
exp(mean(outCM$V0bar)) #229.5738
plot(ts(outCM$tau_V0))
plot(density(outCM$tau_V0))
mean(outCM$tau_V0) #271.2707

plot(ts(outCM$sigmabar))
plot(density(outCM$sigmabar))
exp(mean(outCM$sigmabar)) #0.1542663
plot(ts(outCM$tau_sigma))
plot(density(outCM$tau_sigma))
mean(outCM$tau_sigma) #926.1785

plot(ts(exp(outCM$mat[,,1])))

sigtuneCM = array(0, dim = c(3,3,dim(dfdataCM)[2]))

for(i in 1:dim(dfdataCM)[2]){
sigtuneCM[,,i]=round((1/3)*(2.38^2)*var(outCM$mat[2001:10000,,i]),7)
}

outCM=mhWreMult2(iters=100000, tune = sigtuneCM, y=dfdataCM, times = dftimeCM, init=c(0.31,250.4,1), hypinit = c(0.28, 5.37, 230, 270, 0.15, 930), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Creating matrix of NSCLC tumor data

nmuNSCLC = nmu[nmu[,1]=="NSCLC",]

l1NSCLC=which(nmuNSCLC[,4]==0)
l2NSCLC=which(nmuNSCLC[,4]==0)[-1]

length(l1NSCLC) #29
l1NSCLC

ldiffNSCLC=l2NSCLC-l1NSCLC[-length(l1NSCLC)]
max(ldiffNSCLC) #27

dfdataNSCLC = matrix(NA, nrow=27,ncol = 29)
dftimeNSCLC = matrix(NA, nrow=27,ncol = 29)

for(i in 1:28){
dfdataNSCLC[1:ldiffNSCLC[i],i] = as.matrix(nmuNSCLC[l1NSCLC[i]:(l1NSCLC[(i+1)]-1),3])
dftimeNSCLC[1:ldiffNSCLC[i],i] = as.matrix(nmuNSCLC[l1NSCLC[i]:(l1NSCLC[(i+1)]-1),4])
}

dfdataNSCLC[1:3,29] = as.matrix(nmuNSCLC[286:288,3])
dftimeNSCLC[1:3,29] = as.matrix(nmuNSCLC[286:288,4])

#random effects on NSCLC data

#NB tumour 13 appears to die off and so does not fit this model, hence shall be removed

dfdataNSCLC = dfdataNSCLC[,-13]
dftimeNSCLC = dftimeNSCLC[,-13]

outNSCLC=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,28)), y=dfdataNSCLC, times = dftimeNSCLC, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outNSCLC$alphabar))
plot(density(outNSCLC$alphabar))
exp(mean(outNSCLC$alphabar)) #0.2106894
plot(ts(outNSCLC$tau_alpha))
plot(density(outNSCLC$tau_alpha))
mean(outNSCLC$tau_alpha) #1.994435

plot(ts(outNSCLC$V0bar))
plot(density(outNSCLC$V0bar))
exp(mean(outNSCLC$V0bar)) #219.7177
plot(ts(outNSCLC$tau_V0))
plot(density(outNSCLC$tau_V0))
mean(outNSCLC$tau_V0) #133.5121

plot(ts(outNSCLC$sigmabar))
plot(density(outNSCLC$sigmabar))
exp(mean(outNSCLC$sigmabar)) #0.179478
plot(ts(outNSCLC$tau_sigma))
plot(density(outNSCLC$tau_sigma))
mean(outNSCLC$tau_sigma) #893.5831

plot(ts(exp(outNSCLC$mat[,,1])))

sigtuneNSCLC = array(0, dim = c(3,3,dim(dfdataNSCLC)[2]))

for(i in 1:dim(dfdataNSCLC)[2]){
sigtuneNSCLC[,,i]=round((1/3)*(2.38^2)*var(outNSCLC$mat[2001:10000,,i]),7)
}

outNSCLC=mhWreMult2(iters=100000, tune = sigtuneNSCLC, y=dfdataNSCLC, times = dftimeNSCLC, init=c(0.31,250.4,1), hypinit = c(0.21, 2.00, 220, 130, 0.18, 890), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Creating matrix of PDAC tumor data

nmuPDAC = nmu[nmu[,1]=="PDAC",]

l1PDAC=which(nmuPDAC[,4]==0)
l2PDAC=which(nmuPDAC[,4]==0)[-1]

length(l1PDAC) #38
l1PDAC

ldiffPDAC=l2PDAC-l1PDAC[-length(l1PDAC)]
max(ldiffPDAC) #27

dfdataPDAC = matrix(NA, nrow=27,ncol = 38)
dftimePDAC = matrix(NA, nrow=27,ncol = 38)

for(i in 1:37){
dfdataPDAC[1:ldiffPDAC[i],i] = as.matrix(nmuPDAC[l1PDAC[i]:(l1PDAC[(i+1)]-1),3])
dftimePDAC[1:ldiffPDAC[i],i] = as.matrix(nmuPDAC[l1PDAC[i]:(l1PDAC[(i+1)]-1),4])
}

dfdataPDAC[1:9,38] = as.matrix(nmuPDAC[472:480,3])
dftimePDAC[1:9,38] = as.matrix(nmuPDAC[472:480,4])

#random effects on PDAC data

outPDAC=mhWreMult2(iters=10000, tune=array(diag(c(0.01,0.01,0.01)), dim =c(3,3,38)), y=dfdataPDAC, times = dftimePDAC, init=c(0.31,250.4,1), hypinit = rep(0.01,6), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

plot(ts(outPDAC$alphabar))
plot(density(outPDAC$alphabar))
exp(mean(outPDAC$alphabar)) #0.1513842
plot(ts(outPDAC$tau_alpha))
plot(density(outPDAC$tau_alpha))
mean(outPDAC$tau_alpha) #2.697554

plot(ts(outPDAC$V0bar))
plot(density(outPDAC$V0bar))
exp(mean(outPDAC$V0bar)) #276.5211
plot(ts(outPDAC$tau_V0))
plot(density(outPDAC$tau_V0))
mean(outPDAC$tau_V0) #59.09194

plot(ts(outPDAC$sigmabar))
plot(density(outPDAC$sigmabar))
exp(mean(outPDAC$sigmabar)) #0.1283313
plot(ts(outPDAC$tau_sigma))
plot(density(outPDAC$tau_sigma))
mean(outPDAC$tau_sigma) #429.3068

plot(ts(exp(outPDAC$mat[,,1])))

sigtunePDAC = array(0, dim = c(3,3,dim(dfdataPDAC)[2]))

for(i in 1:dim(dfdataPDAC)[2]){
sigtunePDAC[,,i]=round((1/3)*(2.38^2)*var(outPDAC$mat[2001:10000,,i]),7)
}

outPDAC=mhWreMult2(iters=100000, tune = sigtunePDAC, y=dfdataPDAC, times = dftimePDAC, init=c(0.31,250.4,1), hypinit = c(0.15, 2.7, 280, 59, 0.13, 430), me_a = log(0.2), d_a = 0.0001, g_a = 10, h_a = 0.1, me_s = 0, d_s = 0.0001, g_s = 10, h_s = 0.01, me_V = log(100), d_V = 0.0001, g_V = 10, h_V = 0.01)

#Comparing mean growth rates among groups

plot(density(outGC$alphabar), xlim=c(min(outGC$alphabar), max(outNSCLC$alphabar)), ylim = c(0, 5), main = expression("Posterior densities of the " ~ bar(alpha)), xlab=expression(bar(alpha)))
lines(density(outCRC$alphabar), col = 2)
lines(density(outBRCA$alphabar), col = 3)
lines(density(outPDAC$alphabar), col = 4)
lines(density(outNSCLC$alphabar), col = 5)
lines(density(outCM$alphabar), col = 6)

legend("topleft", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

#Comparing V0bar

plot(density(outGC$V0bar), xlim=c(min(outNSCLC$V0bar), max(outGC$V0bar)), ylim = c(0, 16), main = expression("Posterior densities of the " ~ bar(V[0])), xlab=expression(bar(V[0])))
lines(density(outCRC$V0bar), col = 2)
lines(density(outBRCA$V0bar), col = 3)
lines(density(outPDAC$V0bar), col = 4)
lines(density(outNSCLC$V0bar), col = 5)
lines(density(outCM$V0bar), col = 6)

legend("topleft", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

#Comparing sigmabar

plot(density(outGC$sigmabar), xlim=c(min(outBRCA$sigmabar), -1), ylim = c(0, 12), main = expression("Posterior densities of the " ~ bar(sigma)), xlab=expression(bar(sigma)))
lines(density(outCRC$sigmabar), col = 2)
lines(density(outBRCA$sigmabar), col = 3)
lines(density(outPDAC$sigmabar), col = 4)
lines(density(outNSCLC$sigmabar), col = 5)
lines(density(outCM$sigmabar), col = 6)

legend("topright", legend = c("GC", "CRC", "BRCA", "PDAC", "NSCLC", "CM"), col=seq(1:6), lty=1, cex=0.75)

#Plotting individual solutions
plot(dftimeGC[,1], dfdataGC[,1], type = "n", main = "GC tumour", xlab = "Time", ylab = "Volume")
for(i in 1:1000){
  odeparsW = c(alpha = exp(outGC$mat[i*100,1,1]), rd = ((3/(4*pi))*exp(outGC$mat[i*100,2,1]))^(1/3))
  odeinitW <- c(x = exp(outGC$mat)[i*100,2,1])
  odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
  lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
points(dftimeGC[,1], dfdataGC[,1])

#Bayesian means
mean(exp(outGC$mat[,1,1])) #0.2119637
mean(exp(outGC$mat[,2,1])) #217.865
odeparsW = c(alpha = 0.2119637, rd = ((3/(4*pi))*217.865)^(1/3))
odeinitW <- c(x = 217.865)
odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = 2, lty=2)

#Frequentist means
odeparsW = c(alpha = 3.051384362, rd = 52.403884)
odeinitW <- c(x = 209.4769908909)
odeOutW = ode(odeinitW, na.omit(dftimeGC[,1]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeGC[,1]), odeOutW[,2], col = 4, lty=3)

legend("topleft", legend = c("Bayesian estimate", "Frequentist estimate"), col=c(2,4), lty=c(2,3), cex=0.75)


#attempting to write the MCMC output to data files
write.csv(outGC$alphabar,file="../data/gcAlphaBar.csv")
write.csv(outCRC$alphabar,file="../data/crcAlphaBar.csv")
write.csv(outBRCA$alphabar,file="../data/brcaAlphaBar.csv")
write.csv(outPDAC$alphabar,file="../data/pdacAlphaBar.csv")
write.csv(outNSCLC$alphabar,file="../data/nsclcAlphaBar.csv")
write.csv(outCM$alphabar,file="../data/cmAlphaBar.csv")

write.csv(outGC$V0bar,file="../data/gcV0Bar.csv")
write.csv(outCRC$V0bar,file="../data/crcV0Bar.csv")
write.csv(outBRCA$V0bar,file="../data/brcaV0Bar.csv")
write.csv(outPDAC$V0bar,file="../data/pdacV0Bar.csv")
write.csv(outNSCLC$V0bar,file="../data/nsclcV0Bar.csv")
write.csv(outCM$V0bar,file="../data/cmV0Bar.csv")

write.csv(outGC$sigmabar,file="../data/gcSigmaBar.csv")
write.csv(outCRC$sigmabar,file="../data/crcSigmaBar.csv")
write.csv(outBRCA$sigmabar,file="../data/brcaSigmaBar.csv")
write.csv(outPDAC$sigmabar,file="../data/pdacSigmaBar.csv")
write.csv(outNSCLC$sigmabar,file="../data/nsclcSigmaBar.csv")
write.csv(outCM$sigmabar,file="../data/cmSigmaBar.csv")

write.csv(outGC$mat[,,1],file="../data/gcOut1.csv")

########### From this point we have code that was attempted but ultimately didn't go anywhere ###########

#List of GC and CRC matrices

GC_CRC_data = list(dfdataGC, dfdataCRC)
GC_CRC_time = list(dftimeGC, dftimeCRC)


#THE BIG TEST: fixed and random effects - has been concluded this is not the best way to go

fixtest = mhFixRan(init = c(0.25,250,0.2), y = GC_CRC_data, times = GC_CRC_time, g_a = 10, h_a = 0.1, me_s = log(0.2), g_s = 10, h_s = 0.1, me_V = log(250), g_V = 10, h_V = 0.1)

plot(ts(fixtest$alphabar))
plot(density(fixtest$alphabar))
plot(ts(fixtest$tau_alpha))
plot(density(test2$tau_alpha))

plot(ts(fixtest$B[1,]))
plot(density(fixtest$B[1,]))
var(fixtest$B[1,]) #0.002051034

plot(ts(fixtest$V0bar))
plot(density(fixtest$V0bar))
plot(ts(fixtest$tau_V0))
plot(density(fixtest$tau_V0))

plot(ts(fixtest$sigmabar))
plot(density(fixtest$sigmabar))
plot(ts(fixtest$tau_sigma))
plot(density(fixtest$tau_sigma))

plot(ts(exp(fixtest$mat[,,1])))

sigtunefixtest=round((1/3)*(2.38^2)*var(fixtest$mat[200:1000,,1]),7)

fixtest = mhFixRan(iters = 100000, tune = sigtunefixtest, Btune = 0.002, init = c(0.2,250,0.15), y = GC_CRC_data, times = GC_CRC_time, g_a = 10, h_a = 0.1, me_s = log(0.2), g_s = 10, h_s = 0.1, me_V = log(250), g_V = 10, h_V = 0.1)

plot(ts(fixtest$alphabar[2001:100000]))
plot(density(fixtest$alphabar[2001:100000]))
plot(ts(fixtest$tau_alpha[2001:100000]))
plot(density(test2$tau_alpha[2001:100000]))

plot(ts(fixtest$B[1,]))
plot(density(fixtest$B[1,]))
var(fixtest$B[1,]) #0.0001748487

plot(ts(fixtest$V0bar[2001:100000]))
plot(density(fixtest$V0bar[2001:100000]))
plot(ts(fixtest$tau_V0[2001:100000]))
plot(density(fixtest$tau_V0[2001:100000]))

plot(ts(fixtest$sigmabar[2001:100000]))
plot(density(fixtest$sigmabar[2001:100000]))
plot(ts(fixtest$tau_sigma[2001:100000]))
plot(density(fixtest$tau_sigma[2001:100000]))

plot(ts(exp(fixtest$mat[,,1])))

#plotting individual solution trajectories
plot(dftimeCRC[,45], dfdataCRC[,45], type = "n", main = "Random and fixed effects")
for(i in 1:10000){
odeparsW = c(alpha = exp(fixtest$mat[i*10,1,90]), rd = ((3/(4*pi))*exp(fixtest$mat[i*10,2,90]))^(1/3))
odeinitW <- c(x = exp(fixtest$mat[i*10,2,90]))
odeOutW = ode(odeinitW, na.omit(dftimeCRC[,45]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeCRC[,45]), odeOutW[,2], col = gray(0.8, alpha = 0.2))
}
lines(dftimeCRC[,45], dfdataCRC[,45], type = "p") #Actually looks ok

#Overlaying frequentist estimates where appropriate
odeparsInd = c(alpha = 0.08847871778, rd = 0.00064929)
odeinitInd <- c(x = 180.147878)
odeOutInd = ode(odeinitInd, na.omit(dftimeCRC[,45]), odeFuncW, odeparsW, rtol=1e-3, atol=1e-3)
lines(na.omit(dftimeCRC[,45]), odeOutInd[,2], col = 2, lty = 2)


#Try again with smaller B tuning

fixtest = mhFixRan(iters = 100000, tune = sigtunefixtest, Btune = 0.0002, init = c(0.2,250,0.15), y = GC_CRC_data, times = GC_CRC_time, g_a = 10, h_a = 0.1, me_s = log(0.2), g_s = 10, h_s = 0.1, me_V = log(250), g_V = 10, h_V = 0.1)
