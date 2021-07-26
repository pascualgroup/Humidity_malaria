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

x <- as.numeric(Sys.getenv("PBS_ARRAYID"))
seqP <- seq(from = 1, to = 50000, by = 100)
now.num <- 1#seqP[x]

# here we run the grid parameter function to generate the parameter grid

source("/Users/oscarmauriciosantos/Documents/GitHub/Humidity_malaria/fitting/param_creatGrid.R")

y <- s1 #read.csv("~/Dropbox/paper_humidity_delay/paramGrid_ahmedabad.csv")
param <- as.numeric(y[now.num,])
param.names <- colnames(y)
names(param) <- param.names

source("/Users/oscarmauriciosantos/Documents/GitHub/Humidity_malaria/fitting/pomp_object_RH.R")

# cerate the random walk verctor

rdd<-rw.sd(sigOBS=0.03,
              sigPRO=0.03,
              muS2S1=0.03,
              muEI1=0.03,
              muI2S2=0.03,
              muI2S2=0.03,
              delta=0.03,
              rho=ivp(0.03),
              tau=0.03,
              betaOUT=ivp(0.00001),
              S1_0=ivp(0.03),
              E_0=ivp(0.03),
              I1_0=ivp(0.03),
              S2_0=ivp(0.03),
              I2_0=ivp(0.03),
              K_0=ivp(0.03),
              F_0=ivp(0.03),
              b1=0.03,
              b2=0.03,
              b3=0.03,
              b4=0.03,
              b5=0.03,
              b6=0.03,
              bH=0.03
)


for (i in 1:10000){

seed <- ceiling(runif(1,min=1,max=2^30))
set.seed(seed)
param <- as.numeric(y[i,])
names(param) <- colnames(y)

    tryCatch(mif2(po,
                               Np=1000,
                               Nmif=50,
                               cooling.type="geometric",
                               cooling.fraction.50=0.5,
                               transform=TRUE,
                               start=param,
                               rw.sd=rdd,
                               pred.mean=TRUE,
                               filter.mean=TRUE,max.fail=500),  error = function(e) e) -> mifout
    
    if(length(coef(mifout)) > 0){
      loglik.mif <- replicate(n=5,logLik(pfilter(po,
                                                 params=coef(mifout),Np=1000,max.fail=500)))
      bl <- logmeanexp(loglik.mif,se=TRUE)
      loglik.mif.est <- bl[1]
      loglik.mif.se <- bl[2]
      cat(i,loglik.mif.est,"\n")
      
      ### SAVE OUTPUT
      if(is.finite(loglik.mif.est)){
        par.out <- coef(mifout)
        names(par.out) <- param.names
        if(file.exists("mifOutput_RH.csv")){
          write.table(t(as.matrix(c(i,seed,par.out,loglik.mif.est,loglik.mif.se))), 
                      "mifOutput_RH.csv", append = T, 
                      col.names=F, row.names=F, sep = ",")
        } else{ 
          write.table(t(as.matrix(c(i,seed,par.out,loglik.mif.est,loglik.mif.se))), 
                      "mifOutput_RH.csv", append = T, col.names=c("run","seed", 
                                                                      param.names, "loglik", "loglik.se"), row.names=F, sep = ",") }
      }
    }
  }


