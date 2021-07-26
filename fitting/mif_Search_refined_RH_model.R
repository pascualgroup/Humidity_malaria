#####################
rm(list=ls())
require(pomp)
require(plyr)

dat1 <- read.csv("~/Dropbox/paper_humidity_delay/paper/figure_4/ahmedabad/data_mal.csv",header=T)       # read in data file with the malaria data that is on the data folder
dat2 <- read.csv("~/Dropbox/paper_humidity_delay/paper/figure_4/ahmedabad/covariate3.csv", header=T) #read in data file with covariates that is on the data folder
dat <- join(dat1,dat2)

## data rearreangement and creation of the interannual covariates
splBoth <- loess(dat$pop ~ dat$time)$fitted #
splBoth <- smooth.spline(x=dat$time, splBoth)

var<-cbind(rowMeans(matrix( rep( t( matrix(dat2$RH2,ncol=12,byrow = T)[,3:8] ) , 12) , ncol =  ncol(matrix(dat2$RH2,ncol=12,byrow = T)[,3:8]) , byrow = TRUE )),rep(1:18,12))
dat<-mutate(dat, RH_5=var[ order(var[,2]), 1] )

# here we created a data set 
dat <- mutate(dat, temp = (temp-min(temp))/(max(temp) - min(temp)),
              RH = (RH_5-min(RH_5))/(max(RH_5) - min(RH_5)),
              pop = predict(splBoth, deriv=0, x=dat$time)$y, 
              dpopdt = predict(splBoth, deriv=1, x=dat$time)$y)

now.num <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

y <- read.csv("mifOutput_RH.csv") 
y <- subset(y, !is.na(loglik))
y <- arrange(y, -loglik)
param <- as.numeric(y[now.num,(1:26)])
param.names <- colnames(y)[1:26]
names(param) <- param.names

source("pomp_object_RH.R")


# cerate the random walk verctor

rdd<-rw.sd(sigOBS=0.001,
           sigPRO=0.001,
           muS2S1=0.001,
           muEI1=0.001,
           muI2S2=0.001,
           muI2S2=0.001,
           delta=0.001,
           rho=ivp(0.001),
           tau=0.001,
           betaOUT=ivp(0.00001),
           S1_0=ivp(0.001),
           E_0=ivp(0.001),
           I1_0=ivp(0.001),
           S2_0=ivp(0.001),
           I2_0=ivp(0.001),
           K_0=ivp(0.001),
           F_0=ivp(0.001),
           b1=0.001,
           b2=0.001,
           b3=0.001,
           b4=0.001,
           b5=0.001,
           b6=0.001,
           bH=0.001
)


#########################
for(i in 1:100){
  seed <- ceiling(runif(1,min=1,max=2^30))
  set.seed(seed)
  
  tryCatch(mif2(po,
                Np=5000,
                Nmif=50,
                cooling.type="geometric",
                cooling.fraction.50=0.7,
                transform=TRUE,
                start=param,
                rw.sd=rdd,
                pred.mean=TRUE,
                filter.mean=TRUE,max.fail=500),  error = function(e) e) -> mifout
  mifout <- continue(mifout, Nmif=50, cooling.fraction=0.5)
  mifout <- continue(mifout, Nmif=50, cooling.fraction=0.3)
  
  if(length(coef(mifout)) > 0){
    loglik.mif <- replicate(n=10, logLik(pfilter(po, params=coef(mifout), Np=5000, max.fail=500)))
    bl <- logmeanexp(loglik.mif,se=TRUE)
    loglik.mif.est <- bl[1]
    loglik.mif.se <- bl[2]
    
    if(loglik.mif.est > -1700){
      par.out <- coef(mifout)
      names(par.out) <- param.names
      if(file.exists("mifOutput_RH_refined.csv")){
        write.table(t(as.matrix(c(now.num, seed, par.out, loglik.mif.est, loglik.mif.se))),
                    "mifOutput_RH_refined.csv", append = T, col.names=F, row.names=F, sep = ",")
      } else{
        write.table(t(as.matrix(c(now.num, seed, par.out, loglik.mif.est, loglik.mif.se))),
                    "mifOutput_RH_refined.csv", append = T,
                    col.names=c("run", "seed", param.names, "loglik", "loglik.se"), row.names=F, sep = ",")
      }
    }
  }
}