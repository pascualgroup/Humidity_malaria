########################################
#### Example of a profile iteration ####
########################################

require(pomp)
require(plyr)

dat <- read.csv("DataCovar_1975-2000.csv",header=T)	# read in data file with the malaria data and covariates
y <- read.csv("mifOutputProf_rho.csv", header=T)
y <- arrange(y, -loglik)
y <- y[,2:27]

param <- as.numeric(y[1,])
param.names <- colnames(y)
names(param) <- param.names
source("poObject_TmeanB.R")

index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
temp <- seq(from = 0.0025, to=0.0065,length.out = 201)
param['rho'] <- temp[index]
tempfilt <- -10000

for(i in 1:4){
  seed <- ceiling(runif(1,min=1,max=2^30))
  set.seed(seed)
  
  tryCatch(mif2(
    po,
    Nmif = 150,
    start = param,
    Np = 3000,
    cooling.type="hyperbolic",
    cooling.fraction.50 = 0.3,
    rw.sd = rw.sd(muIS = 0.02, muIQ = 0.02, muQS = 0.02, sigOBS = 0.02, sigPRO = 0.02, muEI = 0.02,
                  tau = 0.02, betaOUT = 0.02, bT4 = 0.02, bT6 = 0.02, b1 = 0.02, b2 = 0.02, 
                  b3 = 0.02, b4 = 0.02, b5 = 0.02, b6 = 0.02, q0 = 0.03,
                  S_0 = 0.03, E_0 = 0.03, I_0 = 0.03, Q_0 = 0.03, K_0 = 0.03, F_0 = 0.03),
    transform=TRUE),  error = function(e) e) -> mifout
  
  if(length(coef(mifout)) > 0){
    loglik.mif <- replicate(n=10, logLik(pfilter(po, params=coef(mifout), Np=3000, max.fail=500)))
    bl <- logmeanexp(loglik.mif,se=TRUE)
    
    if(bl[1] > tempfilt){
      par.out <- coef(mifout)
      names(par.out) <- param.names
      loglik.mif.est <- bl[1]
      loglik.mif.se <- bl[2]
      tempfilt <- loglik.mif.est
      finalSeed <- seed
    }
  }
}

if(file.exists("mifOutputProf_rho.csv")) {
  write.table(t(as.matrix(c(index,par.out,loglik.mif.est,loglik.mif.se))), "mifOutputProf_rho.csv", 
              append = T, col.names=F, row.names=F, sep = ",")
} else{
  write.table(t(as.matrix(c(index,par.out,loglik.mif.est,loglik.mif.se))), "mifOutputProf_rho.csv", 
              append = T, col.names=c("run",param.names, "loglik", "loglik.se"), row.names=F, sep = ",")
}

###############################