#####################
rm(list=ls())
require(pomp)
require(plyr)

dat <- read.csv("DataCovar_1975-2000.csv",header=T)	# read in data file with the malaria data and covariates
now.num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

y <- read.csv("mifOutput_TmeanB.csv") 
y <- subset(y, !is.na(loglik))
y <- arrange(y, -loglik)
param <- as.numeric(y[now.num,(1:26)])
param.names <- colnames(y)[1:26]
names(param) <- param.names
source("poObject_TmeanB.R")

#########################
for(i in 1:3){
  seed <- ceiling(runif(1,min=1,max=2^30))
  set.seed(seed)
  
  tryCatch(mif2(
    po,
    Nmif = 50,
    start = param,
    Np = 3000,
    cooling.type="hyperbolic",
    cooling.fraction.50 = 0.7,
    rw.sd = rw.sd(muEI = 0.03, muIQ = 0.03, muIS = 0.03, muQS = 0.03,  sigPRO = 0.03, sigOBS = 0.03,
                  tau = 0.03, rho = 0.03, betaOUT = 0.03, bT4 = 0.03, bT6 = 0.03, b1 = 0.03, b2 = 0.03, 
                  b3 = 0.03, b4 = 0.03, b5 = 0.03, b6 = 0.03, q0 = 0.03,
                  S_0 = 0.03, E_0 = 0.03, I_0 = 0.03, Q_0 = 0.03 , K_0 = 0.03, F_0 = 0.03),
    transform=TRUE),  error = function(e) e) -> mifout
  mifout <- continue(mifout, Nmif=50, cooling.fraction=0.5)
  mifout <- continue(mifout, Nmif=50, cooling.fraction=0.3)
  
  if(length(coef(mifout)) > 0){
    loglik.mif <- replicate(n=10, logLik(pfilter(po, params=coef(mifout), Np=3000, max.fail=500)))
    bl <- logmeanexp(loglik.mif,se=TRUE)
    loglik.mif.est <- bl[1]
    loglik.mif.se <- bl[2]
    
    if(loglik.mif.est > -1700){
      par.out <- coef(mifout)
      names(par.out) <- param.names
      if(file.exists("mifOutputRef1_TmeanB.csv")){
        write.table(t(as.matrix(c(now.num, seed, par.out, loglik.mif.est, loglik.mif.se))),
                    "mifOutputRef1_TmeanB.csv", append = T, col.names=F, row.names=F, sep = ",")
      } else{
        write.table(t(as.matrix(c(now.num, seed, par.out, loglik.mif.est, loglik.mif.se))),
                    "mifOutputRef1_TmeanB.csv", append = T,
                    col.names=c("run", "seed", param.names, "loglik", "loglik.se"), row.names=F, sep = ",")
      }
    }
  }
}