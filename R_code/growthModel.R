library(deSolve)

#PK-PD model

odeFunc <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dqp <- -(Kpe+Kpt)*qp+Ktp*qT
    dqT <- Kpt*qp-Ktp*qT
    dcp <- (-(Kpe+Kpt)*qp+Ktp*qT)/vp

    return(list(c(dqp, dqT, dcp)))
  })
}

param=c(0.382,0.523,0.196,1.3)
odepars = c(Kpe = param[1], Kpt = param[2], Ktp = param[3], vp = param[4])
odeinit <- c(qp = 14, qT = 0, cp = 11)


dt=0.01
odeOut = ode(odeinit, seq(0,24,length=(24/dt+1)), odeFunc, odepars, rtol=1e-3, atol=1e-3)

plot(ts(odeOut[,2:4],start=0,deltat=dt))
#Compare cp plot with Fig 2 in Evans et al. (2014)


#Tumour growth model

odeFunc2 <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    
    dVt <- kG*Vp+kP*Vh-kDN*Vn
    dVn <- kN*Vh-kDN*Vn
    dTs <- kA*Vh
    if(((3*Vt/(4*pi))^(1/3))<Ts)
    {
     dVc <- 0
    }else{
     dVc <- 4*pi*((((3*Vt)/(4*pi))^(1/3)-Ts)^2)*(((1/3)*(3*Vt/(4*pi))^(-2/3))*dVt - dTs)
    }
    if(Vc<Vn){
     dVh <- 0 
    }else{
     dVh <- dVc - dVn
    }

    dVp <- dVt - dVh - dVn

    return(list(c(dVt, dVn, dTs, dVc, dVh, dVp)))
  })
}

param2=c(0.00252, 0.0011, 5270, 21100, 0)
odepars2 = c(kG = param2[1], kP = param2[2], kDN = param2[3], kN = param2[4], kA = param2[5])
odeinit2 <- c(Vt = 0.322, Vn = 0.0335, Ts = 0.628, Vc = 0.0, Vh = 0.0, Vp = 0.2885)


dt=0.01
odeOut2 = ode(odeinit2, seq(0,1224,length=(1224/dt+1)), odeFunc2, odepars2, rtol=1e-6, atol=1e-6)

plot(ts(round(odeOut2[,2:7],5),start=0,deltat=dt))
#Compare first trace with Fig 3 in Evans et al. (2014)


#Euler solve -just to check!

eul=function(T=1224,dt=0.0001,thin=100,param=c(0.00252, 0.0011, 5270, 21100, 0),init=c(0.322,0.0335,0.628,0.0,0.0,0.2885))
{
 d=length(init)
 n=(T/dt)+1
 n2=(T/(thin*dt))+1
 mat=matrix(0,ncol=d,nrow=n2)
 mat[1,]=init
 kG = param[1]; kP = param[2]; kDN = param[3]; kN = param[4]; kA = param[5]
 Vt = init[1]; Vn = init[2]; Ts = init[3]; Vc = init[4]; Vh = init[5]; Vp = init[6]
 for(i in 1:(n-1)){
    Vt <- Vt+(kG*Vp+kP*Vh-kDN*Vn)*dt
    Vn <- Vn+(kN*Vh-kDN*Vn)*dt
    Ts <- Ts+(kA*Vh)*dt
    if(((3*Vt/(4*pi))^(1/3))<Ts)
    {
     Vc <- 0
    }else{
     Vc <- (4/3)*pi*(((3*Vt)/(4*pi))^(1/3) - Ts)^3
    }
    if(Vc<Vn){
     Vh <- 0 
    }else{
     Vh <- Vc - Vn
    }
    Vp <- Vt - Vh - Vn
    if((i%%thin)==0)
    {
     mat[1+(i/thin),]=c(Vt,Vn,Ts,Vc,Vh,Vp)
    }
 }
 mat
}

mat=eul()
plot(ts(round(mat,5),start=0,deltat=0.01))


##Full growth model with cytotoxicity - ignore for moment! Initial values/param values not quite right...or problem with code

odeFunc3 <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    dqp <- -(Kpe+Kpt)*qp+Ktp*qT
    dqT <- Kpt*qp-Ktp*qT
    dcp <- (-(Kpe+Kpt)*qp+Ktp*qT)/vp
    deta <- cp*IC/(IC+cp)^2
    dVt <- kG*Vp+kP*Vh-kDN*Vn-kDA*Va
    dVn <- kN*Vh-kDN*Vn
    dTs <- kA*Vh
    dVa <- eta*kK*Vp-kDA*Va
    if(((3*Vt/(4*pi))^(1/3))<Ts)
    {
     dVc <- 0
    }else{
     dVc <- 4*pi*((((3*Vt)/(4*pi))^(1/3)-Ts)^2)*(((1/3)*(3*Vt/(4*pi))^(-2/3))*dVt - dTs)
    }
    dVac <- dVc+kDN*(Vc-Vac)+kDA*Vac-kP*Vh
    if(Vc<(Vn+Vac)){
     dVh <- 0 
    }else{
     dVh <- dVc - dVn - dVac
    }

    dVp <- dVt - dVh - dVn - dVa

    return(list(c(dqp, dqT, dcp, deta, dVt, dVn, dTs, dVa, dVc, dVac, dVh, dVp)))
  })
}

param3=c(0.382,0.523,0.196,1.3,0.0000723,0.00252,0.0011,5270,21100,0,0.00165,0.0164)
odepars3 = c(Kpe = param[1], Kpt = param[2], Ktp = param[3], vp = param[4], IC = param[5], kG = param3[6], kP = param3[7], kDN = param3[8], kN = param3[9], kA = param3[10], kDA = param3[11], kK = param3[12])
odeinit3 <- c(qp = 14, qT = 0, cp = 11, eta = 0.99999, Vt = 0.322, Vn = 0.0335, Ts = 0.628, Va = 0.0, Vc = 0.0, Vac = 0.0, Vh = 0.0, Vp = 0.2885)


dt=0.01
odeOut3 = ode(odeinit3, c(0,0.01), odeFunc3, odepars3, rtol=1e-6, atol=1e-6)

plot(ts(round(odeOut3[,2:5],5),start=0,deltat=dt))
plot(ts(round(odeOut3[,6:13],5),start=0,deltat=dt))

Euler solve of the above - seems to "work"

eul2=function(T=1224,dt=0.0001,thin=100,param=c(0.382,0.523,0.196,1.3,0.0000723,0.00252,0.0011,5270,21100,0,0.00165,0.0164),init=c(14,0,11,0.99999,0.322,0.0335,0.628,0.0,0.0,0.0,0.0,0.2885))
{
 d=length(init)
 n=(T/dt)+1
 n2=(T/(thin*dt))+1
 mat=matrix(0,ncol=d,nrow=n2)
 mat[1,]=init
 Kpe = param[1]; Kpt = param[2]; Ktp = param[3]; vp = param[4]; IC = param[5]
 kG = param[6]; kP = param[7]; kDN = param[8]; kN = param[9]; kA = param[10]; kDA = param3[11]; kK = param3[12]
 qp = init[1]; qT = init[2]; cp = init[3]; eta = init[4] 
 Vt = init[5]; Vn = init[6]; Ts = init[7]; Va = init[8]; Vc = init[9]; Vac = init[10]; Vh = init[11]; Vp = init[12]
 for(i in 1:(n-1)){
    qp <- qp+(-(Kpe+Kpt)*qp+Ktp*qT)*dt
    qT <- qT+(Kpt*qp-Ktp*qT)*dt
    cp <- cp+((-(Kpe+Kpt)*qp+Ktp*qT)/vp)*dt
    eta <- cp/(IC+cp)
    Vt <- Vt+(kG*Vp+kP*Vh-kDN*Vn-kDA*Va)*dt
    Vn <- Vn+(kN*Vh-kDN*Vn)*dt
    Ts <- Ts+(kA*Vh)*dt 
    Va <- Va+(eta*kK*Vp-kDA*Va)*dt
    if(((3*Vt/(4*pi))^(1/3))<Ts)
    {
     Vc <- 0
     dVc <- 0
    }else{
     Vc <- (4/3)*pi*(((3*Vt)/(4*pi))^(1/3) - Ts)^3
     dVc <- 4*pi*((((3*Vt)/(4*pi))^(1/3)-Ts)^2)*(((1/3)*(3*Vt/(4*pi))^(-2/3))*(kG*Vp+kP*Vh-kDN*Vn-kDA*Va)-(kA*Vh))
    }
    Vac <- Vac+(dVc+kDN*(Vc-Vac)+kDA*Vac-kP*Vh)*dt
    if(Vc<(Vn+Vac)){
     Vh <- 0 
    }else{
     Vh <- Vc - Vn - Vac
    }
    Vp <- Vt - Vh - Vn - Va
    if((i%%thin)==0)
    {
     mat[1+(i/thin),]=c(qp,qT,cp,eta,Vt,Vn,Ts,Va,Vc,Vac,Vh,Vp)
    }
 }
 mat
}

mat2=eul2()
plot(ts(round(mat2[,1:4],5),start=0,deltat=0.01))
plot(ts(round(mat2[,5:12],5),start=0,deltat=0.01))
