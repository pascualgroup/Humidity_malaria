

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
z <- read.csv("~/Downloads/mifOutputFeb.csv")
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
  data=subset(dat,select=c("time","PF")),
  params=param_pf,
  times="time",
  t0=with(dat,2*time[1]-time[2]),
  covar=covariate_table(time=dat$time,
                        temp=dat$temp,
                        RH=dat$RH,
                        pop=dat$pop,
                        dpopdt=dat$dpopdt,
                        season1=dat$season1,
                        season2=dat$season2,
                        season3=dat$season3,
                        season4=dat$season4,
                        season5=dat$season5,
                        season6=dat$season6,
                        times ="time"),#subset(dat, select=c("time","temp","RH","pop","dpopdt",sprintf("season%d",1:6))),
  rprocess = euler(step.fun = simul_pf, delta.t=1/365),
  rmeasure = rmeas_pf,
  dmeasure = dmeas_pf,
  rinit=initlz_pf,
  parameter_trans(toEst = toEst_pf,fromEst = fromEst_pf),
  accumvars = c("W","cases"),
  statenames = c("cases","S1","E","I1","S2","I2","K","F","W"),
  paramnames = c(par_names,vp_names,sp_names)
) -> pomp_pf_S1S2
