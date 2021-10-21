#mif
require(pomp)
require(magrittr)
require(ggplot2)
require(foreach)
source("pomp_Ahmedabad_pf.R")
pomp_S1EI1I2S2_pf<-create_pomp(tfit = 2009)


now.num <- 3
read.csv("mifOutputFeb.csv")%>% 
  subset(rho > 0.01)%>%
dplyr::arrange(-loglik)->z
z <- z[,1:25]
param_pf <- as.numeric(z[now.num,])
names(param_pf) <- colnames(z)


simulate_pomp(pomp_obj1 = pomp_S1EI1I2S2_pf[[1]],pomp_obj2 = pomp_S1EI1I2S2_pf[[2]],coef_p=param_pf)

rdd_pf<-rw.sd(sigOBS=0.001,
             sigPRO=0.001,
             muS2S1=0,
             muEI1=0,
             muI2S2=0,
             muI2S2=0,
             delta=0,
             rho=ivp(0.001),
             tau=0.00001,
             betaOUT=ivp(0.00001),
             S1_0=ivp(0),
             E_0=ivp(0),
             I1_0=ivp(0),
             S2_0=ivp(0),
             I2_0=ivp(0),
             K_0=ivp(0),
             F_0=ivp(0),
             b1=0.00001,
             b2=0.00001,
             b3=0.00001,
             b4=0.00001,
             b5=0.00001,
             b6=0.00001,
             bH=0.00001
)


pomp_S1EI1I2S2_pf[[1]]%>%mif2(Nmif=1,
                    params=param_pf,
                    rw.sd=rdd_pf,
                    Np=1000,
                    cooling.type="geometric",
                    cooling.fraction.50=0.5)->pomp_S1EI1I2S2_pf_miffed


pomp_S1EI1I2S2_pf_miffed%>%mif2(Nmif=1,
                                rw.sd=rdd_pf,
                                Np=1000,
                                cooling.type="geometric",
                                cooling.fraction.50=0.5)->pomp_S1EI1I2S2_pf_miffed



simulate_pomp(pomp_obj1 =pomp_S1EI1I2S2_pf_miffed,pomp_obj2 = pomp_S1EI1I2S2_pf[[2]],coef_p=coef(pomp_S1EI1I2S2_pf_miffed))

theta=unlist(coef(pomp_S1EI1I2S2_pf_miffed))

theta.t <- partrans(pomp_S1EI1I2S2_pf_miffed,theta,"toEst")
theta.nt <- partrans(pomp_S1EI1I2S2_pf_miffed,theta.t,"fromEst")
estpars <- setdiff(names(theta),c("q0","betaOUT","sigmaPRO","sigmaOBS",
                                  "tau","rho","b1",
                                  "b2","b3",
                                  "b4","b5",
                                  "b6","bH"))

theta.t.hi<- theta.t.lo<- theta.nt
theta.t.lo[estpars] <- theta.nt[estpars]-0.1*theta.nt[estpars]
theta.t.hi[estpars]<- theta.nt[estpars]+0.1*theta.nt[estpars]

runif_design(
  lower= partrans(pomp_S1EI1I2S2_pf_miffed,theta.t.lo,"toEst"),
  upper=partrans(pomp_S1EI1I2S2_pf_miffed,theta.t.hi,"toEst"),
  nseq=20
) %>% as.data.frame() -> pr_1
pr_1 <- as.data.frame(pr_1)
pr_1 <- as.data.frame(t(partrans(pomp_S1EI1I2S2_pf_miffed,t(pr_1),"fromEst")))



foreach (p=iterators::iter(pr_1,"row"),
         .combine=rbind,
         .errorhandling="remove",
         .inorder=FALSE,
         .packages=c("pomp","magrittr","reshape2","plyr"),
         .export="parus"
) %dopar%{
  tic <- Sys.time()
  pomp_S1EI1I2S2_pf_miffed %>% 
    
    mif2(params=unlist(p),
         Nmif=1,
         Np=1000,
         cooling.fraction.50=0.5,
         cooling.type="geometric",
         rw.sd=rdd_pf) %>%
    mif2() -> mf
  
  #print(logLik(mf))
  pf <- replicate(10,pfilter(mf,Np=1000))  ## independent particle filters
  pf$ll <- sapply(pf,logLik)
  ll <- logmeanexp(pf$ll,se=TRUE)
 # nfail <- sapply(pf,getElement,"nfail")  ## number of filtering failures
  #print(ll)
  toc <- Sys.time()
  etime <- toc-tic
  units(etime) <- "hours"
  
  data.frame(as.list(coef(mf)),
             loglik = ll[1],
             loglik.se = ll[2])#,
             #nfail.min = min(nfail),
             #nfail.max = max(nfail))
} -> random_designed_res

saveRDS(random_designed_res,file = "par_exploration")

