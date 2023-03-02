pacman::p_load(DEoptim)
Rosenbrock <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  return(100*(x2-x1*x1)^2+(1-x1)^2)
}
Rosenbrock(c(3.7,5))
lower <- c(-10,-10)
upper <- -lower
set.seed(1234)
DEoptim(Rosenbrock, lower, upper)
DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP=100))


url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/Lab05SetupDRF.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"Lab05SetupDRF.R")
file.edit("Lab05SetupDRF.R")

# Grab out models for Snow and TMWB
# https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R
# becomes: 
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"TMWBFuncs.R")
file.edit("TMWBFuncs.R")
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TISnow.R"
# This will grab the solution for last weeks Lab03 Homework
source(url)
download.file(url,"TISnow.R")
file.edit("TISnow.R")


outTMWB=TMWBmodel(TMWBdf=TMWB)
1-NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred)

################################ TMWB ####################
TMWBoptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- x[6]
  x7 <- x[7]
  outTMWB=TMWBmodel(TMWBdf = TMWB,fcres = x1,Z = x2,SFTmp = x3,bmlt6=x4, bmlt12=x5,Tmlt=x6,Tlag=x7)
  return(1-NSE(Yobs = outTMWB$Qmm,Ysim = outTMWB$Qpred))
}
lower <- c(.01,300,1,.1,.1,1,.1)
upper <- c(.6,3000,6,5,5,6,1)
outDEoptim <- DEoptim(TMWBoptFunc, lower, upper, 
                      DEoptim.control(NP = 80,
                                      itermax = 50, F = 1.2, CR = 0.7))

TMWBmodel



TMWB= BasinData
TMWBoptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- x[6]
  x7 <- x[7]
  x8 <- x[8]
  x9 <- x[9]
  outTMWB=TMWBmodel(TMWBdf = TMWB, fcres = x1, FldCap = x2, WiltPt = x3, Z = x4,
                    SFTmp = x5, bmlt6 = x6, bmlt12 = x7, Tmlt = x8, Tlag = x9)
  return(1-NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred))
}
lower <- c(.01, 0.15, 0.01, 300, 1, 1, 0.1, 1, 0)
upper <- c(0.95, 0.95, 0.5, 1000, 6, 5, 5, 6, 1)
outDEoptim <- DEoptim(TMWBoptFunc, lower, upper,
                      DEoptim.control(NP=200, itermax = 20, F = 1.2, CR =0.7))
detach(TMWBdf)

TMWBinit=TMWBmodel(TMWBdf=TMWB , fcres =  0.162316, FldCap = 0.45, WiltPt = 0.15, Z = 686.909459, 
                   SFTmp = 2.595457, bmlt6 = 1.966454, bmlt12 = 1.426039, Tmlt = 1.07782, Tlag = 0.861917)

plot(TMWBinit$date,TMWBinit$Qmm,type = "l",col="black", xlab = "date", ylab="Q (mm)")
lines(TMWBinit$date,TMWBinit$Qpred,col="blue")
legend("topright", c("Qmm", "Qpred"), col = c("black", "blue"),
       lty = 1:2, cex = 0.8)
title(main= myflowgage$gagename)
NSE(TMWBinit$Qmm,TMWBinit$Qpred)

TMWBopt=TMWBmodel(TMWBdf=TMWB , fcres =  0.011610, FldCap = 0.45, WiltPt = 0.15, Z = 2531.872075, 
                  SFTmp = 1.035587, bmlt6 = 3.253562, bmlt12 = 3.267317, Tmlt = 1.336087, Tlag = 0.828302)

plot(TMWBopt$date,TMWBopt$Qmm,type = "l",col="black", xlab = "date", ylab="Q (mm)")
lines(TMWBopt$date,TMWBopt$Qpred,col="blue")
legend("topright", c("Qmm", "Qpred"), col = c("black", "blue"),
       lty = 1:2, cex = 0.8)
title(main= myflowgage$gagename)
NSE(TMWBopt$Qmm,TMWBopt$Qpred)

############################# CN ###########################
url="https://raw.githubusercontent.com/tyler-26-smith/BSE5304Labs/main/R/CNModel.R"
download.file(url,"CNmodel.R")
file.edit("CNmodel.R")

CNoptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- x[6]
  x7 <- x[7]
  
  
  
  outCN=CNmodel(CNmodeldf=TMWB , CNavg = x1 , IaFrac = x2 , fnc_slope = x3 , fnc_aspect = x4 ,
                func_DAWC = x5 , func_z = x6 , fnc_fcres = x7 )
  return(1-NSE(Yobs = outCN$Qmm, Ysim = outCN$Qpred))
}
lower <- c(50, 0.01, 0.01, 0.01, 0.1,  500 , 0.1)
upper <- c(99, 0.25, 0.25, 0.25, 0.35, 1000 , 0.3)
outDEoptim <- DEoptim(CNoptFunc,lower,upper,
                      DEoptim.control(NP=80, itermax=10, F=1.2, CR=0.7))
detach(CNmodeldf)

CNint=CNmodel(CNmodeldf=TMWB ,CNavg=88.991047,IaFrac = 0.100664, fnc_slope=0.203435    ,
              fnc_aspect=0.213971,func_DAWC=0.308152,func_z=817.803664,fnc_fcres=0.105689 )

plot(CNint$date,CNint$Qmm,type = "l",col="black", xlab = "date", ylab="Q (mm)")
lines(CNint$date,CNint$Qpred,col="blue")
legend("topright", c("Qmm", "Qpred"), col = c("black", "blue"),
       lty = 1:2, cex = 0.8)
title(main= myflowgage$gagename)
NSE(CNint$Qmm,CNint$Qpred)

CNopt=CNmodel(CNmodeldf=TMWB ,CNavg=91.113610,IaFrac = 0.031735, fnc_slope=0.025812    ,
              fnc_aspect=0.070653,func_DAWC=0.342583,func_z=838.651323,fnc_fcres=0.109123 )

plot(CNopt$date,CNopt$Qmm,type = "l",col="black", xlab = "date", ylab="Q (mm)")
lines(CNopt$date,CNopt$Qpred,col="blue")
legend("topright", c("Qmm", "Qpred"), col = c("black", "blue"),
       lty = 1:2, cex = 0.8)
title(main= myflowgage$gagename)
NSE(CNopt$Qmm,CNopt$Qpred)