########################################
#### Example of a profile iteration ####
########################################

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
names(param) <- param.names

source("pomp_object_RH.R")

index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
temp <- seq(from = 0.0025, to=0.0065,length.out = 201)
param['beta'] <- temp[index]
tempfilt <- -10000

for(i in 1:50){
  seed <- ceiling(runif(1,min=1,max=2^30))
  set.seed(seed)
  
  tryCatch(mif2(po,
                Np=3000,
                Nmif=150,
                cooling.type="geometric",
                cooling.fraction.50=0.5,
                transform=TRUE,
                start=param,
                rw.sd=rdd,
                pred.mean=TRUE,
                filter.mean=TRUE,max.fail=500),  error = function(e) e) -> mifout
  
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