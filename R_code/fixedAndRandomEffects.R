#Attempting to implement fixed and random effects

# NB y now a list of matrices - one for each treatment group (tumor type)
# times must now be a list of matrices of the same dimension as y

loglikeWMult2=function(param=c(0.05,0.256,0.1), B,y=mouse1,times=times1,rt=1e-3,at=1e-3)
{
 odepars = c(alpha = exp(log(param[1])+B), rd = ((3/(4*pi))*param[2])^(1/3))
 odeinit = c(x = param[2])
 out2 = ode(odeinit, times, odeFuncW, odepars , rtol=rt, atol=at) 
 x=out2[,2] #latent process 
 return(sum(dnorm(log(y),log(x),param[3],log=TRUE)))
}

mhFixRan=function(iters=1000,tune=diag(c(0.01,0.01,0.01)),init=c(0.01,0.256,0.01), hypinit=c(rep(0.01,6),0),y,times,rt=1e-2,at=1e-3, me_a = 0, d_a = 0.1, g_a = 100, h_a = 1, me_s = log(0.5), d_s = 0.1, g_s = 100, h_s = 1, me_V = log(0.001), d_V = 0.1, g_V = 100, h_V = 1, Bsig = 10, Btune = 1)
{
 ptm = proc.time()
 p=length(init)
 M = sapply(y,dim)[2,]
 N = length(y)
 mat=array(0,dim=c(iters,p,sum(M)))
 mat[1,,]=log(init)
 alphabar=vector('numeric', length = iters)
 tau_alpha=vector('numeric', length = iters)
 V0bar=vector('numeric', length = iters)
 tau_V0=vector('numeric', length = iters)
 sigmabar=vector('numeric', length = iters)
 tau_sig=vector('numeric', length = iters)
 B = matrix(nrow = N-1, ncol = iters)
 
 alphabar[1] = log(hypinit[1])
 tau_alpha[1] = hypinit[2]
 V0bar[1] = log(hypinit[3])
 tau_V0[1] = hypinit[4]
 sigmabar[1] = log(hypinit[5])
 tau_sig[1] = hypinit[6]
 B[1,] = hypinit[7]
 curr=matrix(log(init),nrow = p, ncol = sum(M))
 llikecurr=rep(-1e8, sum(M)) #accept first
 count=rep(1, sum(M))
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
 
   for(j in 1:M[1]){
    can=rmvn(curr[,j],tune)
    #B = 0 for first group
    llikecan=loglikeWMult2(exp(can), B = 0, na.omit(y[[1]][,j]),na.omit(times[[1]][,j]),rt,at)
    laprob=llikecan-llikecurr[j]+logprior(exp(can[1]),alphabar[i],1/sqrt(tau_alpha[i]))+logprior(exp(can[2]),V0bar[i],1/sqrt(tau_V0[i]))+logprior(exp(can[3]),sigmabar[i],1/sqrt(tau_sig[i]))-logprior(exp(curr[2,j]),V0bar[i],1/sqrt(tau_V0[i])) -logprior(exp(curr[1,j]),alphabar[i],1/sqrt(tau_alpha[i]))-logprior(exp(curr[3,j]),sigmabar[i],1/sqrt(tau_sig[i]))
    if(log(runif(1))<laprob)
    {
     curr[,j]=can
     llikecurr[j]=llikecan
     count[j]=count[j]+1
    }
    mat[i,,j]=curr[,j]
   }
   
   for(k in 2:N){
   #propose candidate B value and use current values of alpha (and other parameters) to determine acceptance prob.
   Bllikcurr = sum(llikecurr[M[k-1]:(M[k-1]+M[k])])
   Bcan = rnorm(1, B[k-1, i-1], Btune)
   Bllikcan = 0
   for(j in 1:M[k]){
    Bllikcan=Bllikcan + loglikeWMult2(exp(curr[,(M[k-1]+j)]), B = Bcan ,na.omit(y[[k]][,j]),na.omit(times[[k]][,j]),rt,at)
   }
   Blaprob = Bllikcan - Bllikcurr + dnorm(Bcan, 0, Bsig, log = T) - dnorm(B[k-1, i-1], 0, Bsig, log = T)
   if(log(runif(1))<Blaprob)
    {
   B[k-1, i] = Bcan
   }
   else{
   B[k-1, i] = B[k-1, i-1]
   }
   for(j in 1:M[k]){
    can=rmvn(curr[,(M[k-1]+j)],tune)
   #Fix current value of B and use in inference of a
    llikecan=loglikeWMult2(exp(can), B = B[k-1,i] ,na.omit(y[[k]][,j]),na.omit(times[[k]][,j]),rt,at)
   #Unsure if changing B should alter llikecurr, if so, this is inaccurate
    laprob=llikecan-llikecurr[(M[k-1]+j)]+logprior(exp(can[1]),alphabar[i],1/sqrt(tau_alpha[i]))+logprior(exp(can[2]),V0bar[i],1/sqrt(tau_V0[i]))+logprior(exp(can[3]),sigmabar[i],1/sqrt(tau_sig[i]))-logprior(exp(curr[2,j]),V0bar[i],1/sqrt(tau_V0[i]))-logprior(exp(curr[1,j]),alphabar[i],1/sqrt(tau_alpha[i]))-logprior(exp(curr[3,j]),sigmabar[i],1/sqrt(tau_sig[i]))
    if(log(runif(1))<laprob)
    {
     curr[,(M[k-1]+j)]=can
     llikecurr[(M[k-1]+j)]=llikecan
     count[(M[k-1]+j)]=count[(M[k-1]+j)]+1
    }
    mat[i,,(M[k-1]+j)]=curr[,(M[k-1]+j)]
   }
 }
 }
 print(count/iters)
 print(proc.time()-ptm)
 return(list(alphabar = alphabar, tau_alpha = tau_alpha, V0bar = V0bar, tau_V0 = tau_V0, sigmabar = sigmabar, tau_sigma = tau_sig, mat = mat, B = B))
}
