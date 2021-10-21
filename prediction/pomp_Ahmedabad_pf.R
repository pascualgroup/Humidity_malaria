
############ data #################

create_pomp <-function(tfit=2011){
dat1 <- read.csv("data_Ahme_M.csv",header=T)       # read in data file with the malaria data and the mosquito densityies and breeding spots. 
dat2 <- read.csv("covariates_ahmedabad.csv", header=T) #read in data file with covariates  # THIS IS CREATED with create_covariate_file3.R, and uses the temperature, rain relative humidity. 
names(dat2)[1]<-"time"
dat <- plyr::join(dat1,dat2)%>%
  dplyr::filter(time<2015)
#data before and during 2011
dat%>% dplyr::filter(time < (tfit + 1))->dat_fit
dat2%>% dplyr::filter(time < (tfit + 1))->dat_cov_fit

#data after 2011
dat%>% dplyr::filter(time >= (tfit + 1) & time <2015)->dat_predict
dat2%>% dplyr::filter(time >= (tfit + 1)& time <2015)->covar_predict
#larv <- read.csv("~/Dropbox/paper_humidity_delay/paper/figure_4/ahmedabad/larvdat.csv", header=T)

splBoth <- loess(dat_fit$pop ~ dat_fit$time)$fitted #
splBoth <- smooth.spline(x=dat_fit$time, splBoth)

splBoth_p <- loess(dat_predict$pop ~ dat_predict$time)$fitted #
splBoth_p <- smooth.spline(x=dat_predict$time, splBoth_p)


var_tot<-cbind(rowMeans(matrix(rep( t( matrix(dat2$RH,ncol=12,byrow = T)[,3:8] ) , 12),
                           ncol =  ncol(matrix(dat2$RH,ncol=12,byrow = T)[,3:8]),
                           byrow = TRUE )),
               rep(1:18,12)
)


var<-cbind(rowMeans(matrix(rep( t( matrix(dat_cov_fit$RH,ncol=12,byrow = T)[,3:8] ) , 12),
                           ncol =  ncol(matrix(dat_cov_fit$RH,ncol=12,byrow = T)[,3:8]),
                           byrow = TRUE )),
           rep(1:(tfit+1),12)
           )
dat_fit<-dplyr::mutate(dat_fit, RH_5=var[ order(var[1:dim(dat_fit)[1],2]), 1] )

dat_fit <- dplyr::mutate(dat_fit, temp = (temp-min(temp))/(max(temp) - min(temp)),
              RH = (RH_5-min(RH_5))/(max(RH_5) - min(RH_5)),
              pop = predict(splBoth, deriv=0, x=dat_fit$time)$y, 
              dpopdt = predict(splBoth, deriv=1, x=dat_fit$time)$y)

dat_fit<- dplyr::mutate(dat_fit,month=round((time - floor(time))*12+1,digits=0))
dat_fit<-dplyr::mutate(dat_fit,year=floor(time))
seasdat<-dat_fit[,c("month","season1","season2","season3","season4","season5","season6")]
seasdat<-reshape2::melt(seasdat,id.vars = "month")




######################################################
var2<-cbind(rowMeans(matrix(rep(t(matrix(covar_predict$RH,ncol=12,byrow = T)[,3:8]),12),
                           ncol=  ncol(matrix(covar_predict$RH,ncol=12,byrow = T)[,3:8]),byrow = TRUE )),
           rep(1:(2015-tfit-1),12))

dat_predict<-dplyr::mutate(dat_predict, RH_5=var2[ order(var2[1:dim(dat_predict)[1],2]), 1] )

dat_predict <- dplyr::mutate(dat_predict, temp = (temp-min(temp))/(max(temp) - min(temp)),
                         RH = (RH_5-min(RH_5))/(max(RH_5) - min(RH_5)),
                         pop = predict(splBoth_p, deriv=0, x=dat_predict$time)$y, 
                         dpopdt = predict(splBoth_p, deriv=1, x=dat_predict$time)$y)

dat_predict<- dplyr::mutate(dat_predict,month=round((time - floor(time))*12+1,digits=0))
dat_predict<-dplyr::mutate(dat_predict,year=floor(time))




#splines<-ggplot(seasdat,aes(x=factor(month),y=scales::rescale(value),group=variable,colour=factor(variable)))+
#  geom_smooth(span = 0.98, method = 'loess',se=T) +
#  xlab("Month")+ylab("Value")

#ggsave("~/Dropbox//Disertation/figure_7ch4.pdf", splines, width = 10, height = 5)

scales::rescale()

simul_pf <- Csnippet("
                  // compute transmission rate 
  double betaIN = exp(b1*season1+b2*season2+b3*season3+b4*season4+b5*season5+b6*season6+b4*season4*bH*RH);
                  
                  // gamma white noise
                  double dW = rgammawn(sigPRO,dt);
                  
                  // force of infection
                  double foi = (betaOUT+((I1+q0*I2)/pop)*betaIN)*(dW/dt);
                  
                  double dBS1 = (delta*pop+dpopdt)*dt;
                  double dS2S1 = muS2S1*S2*dt;
                  double dS1E = F*S1*dt;
                  
                  double dEI1 = muEI1*E*dt;
                  double dI1S2 = muI1S2*I1*dt; // ask about the fIR
                  double dS2I2 = F*S2*dt;
                  double dI2S2 = muI2S2*I2*dt;
                  double dS1D = delta*S1*dt;
                  double dED = delta*E*dt;
                  double dI1D = delta*I1*dt;
                  double dS2D = delta*S2*dt;
                  double dI2D = delta*I2*dt;
                  double dKK = ((foi-K)/(tau/2.0))*dt;
                  double dFF = ((K-F)/(tau/2.0))*dt;
                  
                  // compute equations
                  S1 += dBS1 + dS2S1 + dS1E - dS1D;
                  E += dS1E - dEI1 - dED;
                  I1 += dEI1 - dI1S2 - dI1D;
                  S2 += dI1S2 - dS2S1 - dS2I2 +dI2S2 - dS2D;
                  I2 += dS2I2 - dI2S2 - dI2D;
                  K += dKK;
                  F += dFF;
                  cases += rho*dEI1;
                  W += (dW-dt)/sigPRO;
                  ")

############ rmeas #################
rmeas_pf <- Csnippet("
                  double size = 1.0/sigOBS/sigOBS;
                  PF = rnbinom_mu(size,cases);
                  ")

############ dmeas #################
dmeas_pf <- Csnippet("
                  double size = 1.0/sigOBS/sigOBS;
                  lik = dnbinom_mu(PF,size,cases+0.1,1);
                  if (!give_log) lik = exp(lik);
                  ")

############ fromEst #################
fromEst_pf <- Csnippet("
                    TsigOBS = expit(sigOBS);
                    TsigPRO = expit(sigPRO);
                    TmuS2S1 = exp(muS2S1);
                    TmuEI1 = exp(muEI1);
                    TmuI1S2 = exp(muI1S2);
                    TmuI2S2 = exp(muI2S2);
                    TbetaOUT = exp(betaOUT);
                    Trho = expit(rho);
                    Ttau = expit(tau);
                    Tq0 = expit(q0);
                    TS1_0 = expit(S1_0);
                    TE_0 = expit(E_0);
                    TI1_0 = expit(I1_0);
                    TS2_0 = expit(S2_0);
                    TI2_0 = expit(I2_0);
                    TK_0 = expit(K_0);
                    TF_0 = expit(F_0);
                    double sum = TS1_0 + TE_0 + TI1_0 + TS2_0 + TI2_0;
                    TS1_0 /= sum;
                    TE_0 /= sum;
                    TI1_0 /= sum;
                    TS2_0 /= sum;
                    TI2_0 /= sum;
                    ")

############ fromEst #################

toEst_pf <- Csnippet("
                  TsigOBS = logit(sigOBS);
                  TsigPRO = logit(sigPRO);
                  TmuS2S1 = log(muS2S1);
                  TmuEI1 = log(muEI1);
                  TmuI1S2 = log(muI1S2);
                  TmuI2S2 = log(muI2S2);
                  TbetaOUT = log(betaOUT);
                  Trho = logit(rho);
                  Ttau = logit(tau);
                  Tq0 = logit(q0);
                  TS1_0 = logit(S1_0);
                  TE_0 = logit(E_0);
                  TI1_0 = logit(I1_0);
                  TS2_0 = logit(S2_0);
                  TI2_0 = logit(I2_0);
                  TK_0 = logit(K_0);
                  TF_0 = logit(F_0);
                  
                  double sum = S1_0 + E_0 + I1_0 + S2_0 + I2_0;
                  TS1_0 /= log(S1_0/sum);
                  TE_0 /= log(E_0/sum);
                  TI1_0 /= log(I1_0/sum);
                  TS2_0 /= log(S2_0/sum);
                  TI2_0 /= log(I2_0/sum);
                  ")

############ initlz #################

initlz_pf <- Csnippet("
                   double m = pop/(S1_0 + E_0 + I1_0 + S2_0 + I2_0);
                   
                   S1 = nearbyint(m*S1_0);
                   E = nearbyint(m*E_0);
                   I1 = nearbyint(m*I1_0);
                   S2 = nearbyint(m*S2_0);
                   I2 = nearbyint(m*I2_0);
                   
                   K = K_0;
                   F = F_0;
                   cases = 0;
                   W = 0;
                   ")


############ indexCluster #################
now.num <- 3
z <- read.csv("mifOutputFeb.csv")
z <- subset(z, rho > 0.01)
z <- dplyr::arrange(z, -loglik)
z <- z[,1:25]
param_pf <- as.numeric(z[now.num,])
names(param_pf) <- colnames(z)


par_names=c("sigOBS","sigPRO","muS2S1","muEI1","muI1S2","muI2S2","betaOUT","delta","rho","tau","q0");
vp_names <-c("S1_0","E_0","I1_0","I2_0","S2_0","K_0","F_0")
sp_names <-c("b1","b2","b3","b4","b5","b6","bH")
############ pomp Object #################
pomp(
    data=subset(dat_fit,select=c("time","PF")),
  params=param_pf,
  times="time",
  t0=with(dat_fit,2*time[1]-time[2]),
  covar=covariate_table(time=dat_fit$time,
                        temp=dat_fit$temp,
                        RH=dat_fit$RH,
                        pop=dat_fit$pop,
                        dpopdt=dat_fit$dpopdt,
                        season1=dat_fit$season1,
                        season2=dat_fit$season2,
                        season3=dat_fit$season3,
                        season4=dat_fit$season4,
                        season5=dat_fit$season5,
                        season6=dat_fit$season6,
                        times ="time"),#subset(dat, select=c("time","temp","RH","pop","dpopdt",sprintf("season%d",1:6))),
  rprocess = euler(step.fun = simul_pf, delta.t=1/365),
  rmeasure = rmeas_pf,
  dmeasure = dmeas_pf,
  rinit=initlz_pf,
  parameter_trans(toEst = toEst_pf,fromEst = fromEst_pf),
  accumvars = c("W","cases"),
  statenames = c("cases","S1","E","I1","S2","I2","K","F","W"),
  paramnames = c(par_names,vp_names,sp_names)
) -> pomp_S1EI1I2S2_pf

pomp(
  data=subset(dat_predict,select=c("time","PF")),
  
  times="time",
  t0=with(dat_predict,2*time[1]-time[2]),
  covar=covariate_table(time=dat_predict$time,
                        temp=dat_predict$temp,
                        RH=dat_predict$RH,
                        pop=dat_predict$pop,
                        dpopdt=dat_predict$dpopdt,
                        season1=dat_predict$season1,
                        season2=dat_predict$season2,
                        season3=dat_predict$season3,
                        season4=dat_predict$season4,
                        season5=dat_predict$season5,
                        season6=dat_predict$season6,
                        times ="time"),#subset(dat, select=c("time","temp","RH","pop","dpopdt",sprintf("season%d",1:6))),
  rprocess = euler(step.fun = simul_pf, delta.t=1/365),
  rmeasure = rmeas_pf,
  dmeasure = dmeas_pf,
  rinit=initlz_pf,
  parameter_trans(toEst = toEst_pf,fromEst = fromEst_pf),
  accumvars = c("W","cases"),
  statenames = c("cases","S1","E","I1","S2","I2","K","F","W"),
  paramnames = c(par_names,vp_names,sp_names)
) -> pomp_S1EI1I2S2_pf_sim

res=list(pomp_S1EI1I2S2_pf,pomp_S1EI1I2S2_pf_sim)
return(res)
}

###################

simulate_pomp<-function(pomp_obj1,pomp_obj2,coef_p){
 


  
  pfout <- pfilter(pomp_obj1,params=coef_p,Np=1000,save.states=T)
  S1_0 <- median(pfout@saved.states[[length(pfout@saved.states)]]['S1',])
  E0<-median(pfout@saved.states[[length(pfout@saved.states)]]['E',])
  I1_0 <- median(pfout@saved.states[[length(pfout@saved.states)]]['I1',])
  S2_0 <- median(pfout@saved.states[[length(pfout@saved.states)]]['S2',])
  I2_0 <- median(pfout@saved.states[[length(pfout@saved.states)]]['I2',])
  
  coef(pomp_obj1)[19]=S1_0/(S1_0+E0+I1_0+S2_0+I2_0)
  coef(pomp_obj1)[20]=E0/(S1_0+E0+I1_0+S2_0+I2_0)
  coef(pomp_obj1)[21]=I1_0/(S1_0+E0+I1_0+S2_0+I2_0)
  coef(pomp_obj1)[22]=S2_0/(S1_0+E0+I1_0+S2_0+I2_0)
  coef(pomp_obj1)[23]=I2_0/(S1_0+E0+I1_0+S2_0+I2_0)
  coef(pomp_obj1)[24]=median(pfout@saved.states[[length(pfout@saved.states)]]['K',])
  coef(pomp_obj1)[25]=median(pfout@saved.states[[length(pfout@saved.states)]]['F',])
  
  
  pomp_obj2%>%
    simulate(nsim=300,
             params=coef(pomp_obj1),
             include.data=TRUE,
             format= "data.frame") %>%
    subset(time>1997,select=c(time,.id,PF)) %>%
    dplyr::mutate(data=.id=="data") %>%
    plyr::ddply(~time+data,dplyr::summarize,
                p=c(0.25,0.5,0.75),q=quantile(PF,prob=p,names=FALSE,na.rm=T)) %>%
    dplyr::mutate(p=plyr::mapvalues(p,from=c(0.25,0.5,0.75),to=c("lo","cases","hi")),
                  data=plyr::mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
    reshape2::dcast(time+data~p,value.var='q')->sim_obj
#  browser()
    ggplot() + geom_ribbon(data=subset(sim_obj,data=="simulation"),aes(x=time,ymin=lo,ymax=hi),alpha=I(0.2),fill="grey70")+
    geom_line(data=subset(sim_obj,data=="data"),aes(x=time,y=cases),color="red")+
      geom_line(data=subset(sim_obj,data=="simulation"),aes(x=time,y=cases),color="black",linetype=2)+
      theme_bw()+scale_y_continuous("Cases")+
    scale_color_manual("",values=c("red","grey70"))
  
}
