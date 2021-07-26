require(tgp)
s1 <- lhs(10000,rbind(
  c(0,0.2),  # sigOBS
  c(0,0.2),  # sigPRO
  c(0,100), # muS2S1
  c(0,100),  # muEI1
  c(0,100),  # muI1S2
  c(0,100), # muI2S2
  c(0,0.1),    # betaOUT
  c(0,0.1),    # rho
  c(0.2,0.8), #tau
  c(-3,3),   # bH   
  c(0,1),    # q0
  c(-5,5), # b1
  c(-5,5), # b2
  c(-5,5), # b3
  c(-5,5), # b4
  c(-5,5), # b5
  c(-5,5), # b6
  c(0.02,0.02), # delta
  c(0,0.5),  # S.1
  c(0,0.2),  # E.0
  c(0,0.2),  # I.1
  c(0,0.2),  # S.2
  c(0,0.2),    # I.2
  c(0,0.5),     # F.0
  c(0.02,0.5)     # K.0
)
)



colnames(s1) <- c(
  "sigOBS",
  "sigPRO",
  "muS2S1",
  "muEI1",
  "muI1S2",
  "muI2S2",
  "betaOUT",
  "rho",
  "tau",
  "bH",
  "q0",
  "b1",
  "b2",
  "b3",
  "b4",
  "b5",
  "b6",
  "delta",
  "S1_0",
  "E_0",
  "I1_0",
  "S2_0",
  "I2_0",
  "K_0",
  "F_0"
)

s1 <- as.data.frame(s1)

write.csv(s1,file="~/paramGrid_surat.csv",row.names=F)	


