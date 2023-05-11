### TEST performance of the package on simulations

# We will test the package with correctly specificied Q and g functions
library(MargIntTmle)

# rexpit function
rexpit <- function (x) rbinom(length(x), 1, plogis(x))

gen.sim.data <- function(N, do.A1, do.A2) {
  # baseline confounders:
  conf1 <- rbinom(N, size = 1, prob = 0.50)
  conf2 <- rbinom(N, size = 1, prob = 0.20)
  conf3 <- rbinom(N, size = 1, prob = 0.70)

  # exposure
  if (is.null(do.A1) & is.null(do.A2)) {
    fact.A1 <- rexpit(0.5 + log(2) * conf1 + log(0.6) * conf2)
    fact.A2 <- rexpit(-1 + log(2.3) * conf1 + log(1.8) * conf3)
  } else {
    fact.A1 = rep(do.A1, N)
    fact.A2 = rep(do.A2, N)
  }

  # outcome
  Y.bin <- rexpit(-3 + log(0.5) * conf1 + log(2.1) * conf2 + log(0.7) * conf3 +
                    log(1.8) * fact.A1 + log(2.2) * fact.A2 + log(2) * fact.A1 * fact.A2)
  Y.cont <- rnorm(N,
                  mean = 100 + 20 * conf1 - 30 * conf2 + 40 * conf3 +
                    20 * fact.A1 + 40 * fact.A2 + 30 * fact.A1 * fact.A2,
                  sd = 20)

  # data.frame
  data.sim <- data.frame(conf1, conf2, conf3, fact.A1, fact.A2, Y.bin, Y.cont)

  # output: simulated data.frame or mean outcome under counterfactual intervention do(A1,A2)
  if (is.null(do.A1) & is.null(do.A2)) {
    return(data.sim)
  } else {
    return(list(mean.Y.bin = mean(Y.bin),
                mean.Y.cont = mean(Y.cont)))
  }
}

### 0) True values
### 0.1) Estimated by simulations
# set.seed(12345)
# p00 <- gen.sim.data(N = 1e6, do.A1 = 0, do.A2 = 0)
# p10 <- gen.sim.data(N = 1e6, do.A1 = 1, do.A2 = 0)
# p01 <- gen.sim.data(N = 1e6, do.A1 = 0, do.A2 = 1)
# p11 <- gen.sim.data(N = 1e6, do.A1 = 1, do.A2 = 1)
#
# ### For Y binary
# p00$mean.Y.bin # [1] 0.03441
# p10$mean.Y.bin # [1] 0.059572
# p01$mean.Y.bin # [1] 0.071873
# p11$mean.Y.bin # [1] 0.212607
#
# # Risk Difference
# # RD_A1_A2is0
# p10$mean.Y.bin - p00$mean.Y.bin # [1] 0.025162
# # RD_A1_A2is1
# p11$mean.Y.bin - p01$mean.Y.bin # [1] 0.140734
# # RD_A2_A1is0
# p01$mean.Y.bin - p00$mean.Y.bin # [1] 0.037463
# # RD_A2_A1is1
# p11$mean.Y.bin - p10$mean.Y.bin # [1] 0.153035
#
# # Relative risk
# # RR_A1_A2is0
# p10$mean.Y.bin / p00$mean.Y.bin # [1] 1.731241
# # RR_A1_A2is1
# p11$mean.Y.bin / p01$mean.Y.bin # [1] 2.958093
# # RR_A2_A1is0
# p01$mean.Y.bin / p00$mean.Y.bin # [1] 2.088724
# # RR_A2_A1is1
# p11$mean.Y.bin / p10$mean.Y.bin # [1] 3.568908
#
# # interaction effects
# # a.INT
# p11$mean.Y.bin - p10$mean.Y.bin - p01$mean.Y.bin + p00$mean.Y.bin # [1] 0.115572
# # m.INT
# (p11$mean.Y.bin * p00$mean.Y.bin) / (p10$mean.Y.bin * p01$mean.Y.bin)  # [1] 1.708655
# # RERI
# (p11$mean.Y.bin - p10$mean.Y.bin - p01$mean.Y.bin + p00$mean.Y.bin) / p00$mean.Y.bin # [1] 3.358675
#
# ### For Y continuous
# p00$mean.Y.cont # [1] 132.0172
# p10$mean.Y.cont # [1] 152.0218
# p01$mean.Y.cont # [1] 172.0286
# p11$mean.Y.cont # [1] 221.992
#
# # Risk Difference
# # RD_A1_A2is0
# p10$mean.Y.cont - p00$mean.Y.cont # [1] 20.00453
# # RD_A1_A2is1
# p11$mean.Y.cont - p01$mean.Y.cont # [1] 49.96336
# # RD_A2_A1is0
# p01$mean.Y.cont - p00$mean.Y.cont # [1] 40.01139
# # RD_A2_A1is1
# p11$mean.Y.cont - p10$mean.Y.cont # [1] 69.97023
#
# # interaction effects
# # a.INT
# p11$mean.Y.cont - p10$mean.Y.cont - p01$mean.Y.cont + p00$mean.Y.cont # [1] 29.95883


### 0.2) Estimated by analytical computation
S <- cbind(expand.grid(c(0,1),c(0,1),c(0,1)),
           rep(NA,n=2^3),rep(NA,n=2^3),
           rep(NA,n=2^3),rep(NA,n=2^3),
           rep(NA,n=2^3),rep(NA,n=2^3),
           rep(NA,n=2^3),rep(NA,n=2^3))
colnames(S) <- list("conf1","conf2","conf3",
                    "sum.Y00.bin","sum.Y00.cont",
                    "sum.Y10.bin","sum.Y10.cont",
                    "sum.Y01.bin","sum.Y01.cont",
                    "sum.Y11.bin","sum.Y11.cont")
for (k in 1:8) {
  S[k,"sum.Y00.bin"] <- plogis(-3 + log(0.5) * S[k,"conf1"] +
                                 log(2.1) * S[k,"conf2"] +
                                 log(0.7) * S[k,"conf3"] +
                                 log(1.8) * 0 +
                                 log(2.2) * 0 +
                                 log(2) * 0 * 0) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y00.cont"] <- (100 + 20 * S[k,"conf1"] -
                            30 * S[k,"conf2"] +
                            40 * S[k,"conf3"] +
                            20 * 0 +
                            40 * 0 +
                            30 * 0 * 0) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y10.bin"] <- plogis(-3 + log(0.5) * S[k,"conf1"] +
                                 log(2.1) * S[k,"conf2"] +
                                 log(0.7) * S[k,"conf3"] +
                                 log(1.8) * 1 +
                                 log(2.2) * 0 +
                                 log(2) * 1 * 0) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y10.cont"] <- (100 + 20 * S[k,"conf1"] -
                            30 * S[k,"conf2"] +
                            40 * S[k,"conf3"] +
                            20 * 1 +
                            40 * 0 +
                            30 * 1 * 0) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y01.bin"] <- plogis(-3 + log(0.5) * S[k,"conf1"] +
                                 log(2.1) * S[k,"conf2"] +
                                 log(0.7) * S[k,"conf3"] +
                                 log(1.8) * 0 +
                                 log(2.2) * 1 +
                                 log(2) * 0 * 1) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y01.cont"] <- (100 + 20 * S[k,"conf1"] -
                            30 * S[k,"conf2"] +
                            40 * S[k,"conf3"] +
                            20 * 0 +
                            40 * 1 +
                            30 * 0 * 1) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y11.bin"] <- plogis(-3 + log(0.5) * S[k,"conf1"] +
                                 log(2.1) * S[k,"conf2"] +
                                 log(0.7) * S[k,"conf3"] +
                                 log(1.8) * 1 +
                                 log(2.2) * 1 +
                                 log(2) * 1 * 1) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))

  S[k,"sum.Y11.cont"] <- (100 + 20 * S[k,"conf1"] -
                            30 * S[k,"conf2"] +
                            40 * S[k,"conf3"] +
                            20 * 1 +
                            40 * 1 +
                            30 * 1 * 1) *
    (0.50^S[k,"conf1"]) * ((1-0.50)^(1-S[k,"conf1"])) *
    (0.20^S[k,"conf2"]) * ((1-0.20)^(1-S[k,"conf2"])) *
    (0.70^S[k,"conf3"]) * ((1-0.70)^(1-S[k,"conf3"]))
}

## Y binary
# point estimates
p00 <- sum(S[,"sum.Y00.bin"]) # 0.03440588
p10 <- sum(S[,"sum.Y10.bin"]) # 0.05986509
p01 <- sum(S[,"sum.Y01.bin"]) # 0.07198179
p11 <- sum(S[,"sum.Y11.bin"]) # 0.2120109

# Risk Difference
# RD_A1_A2is0
RD.bin_A1_A2is0 <- p10 - p00 # [1] 0.02545921
# RD_A1_A2is1
RD.bin_A1_A2is1 <- p11 - p01 # [1] 0.1400291
# RD_A2_A1is0
RD.bin_A2_A1is0 <- p01 - p00 # [1] 0.03757591
# RD_A2_A1is1
RD.bin_A2_A1is1 <- p11 - p10 # [1] 0.1521458

# Relative risk
# RR_A1_A2is0
RR.bin_A1_A2is0 <- p10 / p00 # [1] 1.739967
# RR_A1_A2is1
RR.bin_A1_A2is1 <- p11 / p01 # [1] 2.94534
# RR_A2_A1is0
RR.bin_A2_A1is0 <- p01 / p00 # [1] 2.092137
# RR_A2_A1is1
RR.bin_A2_A1is1 <- p11 / p10 # [1] 3.541478

# interaction effects
# a.INT
a.INT.bin <- p11 - p10 - p01 + p00 # [1] 0.1145699
# m.INT
m.INT.bin <- (p11 * p00) / (p10 * p01)  # [1] 1.692757
# RERI
RERI.bin <- (p11 - p10 - p01 + p00) / p00 # [1] 3.329951

## Y continuous
# point estimates
E.Y00 <- sum(S[,"sum.Y00.cont"]) # 132
E.Y10 <- sum(S[,"sum.Y10.cont"]) # 152
E.Y01 <- sum(S[,"sum.Y01.cont"]) # 172
E.Y11 <- sum(S[,"sum.Y11.cont"]) # 222

# Risk Difference
# RD_A1_A2is0
RD.cont_A1_A2is0 <- E.Y10 - E.Y00 # [1] 20
# RD_A1_A2is1
RD.cont_A1_A2is1 <- E.Y11 - E.Y01 # [1] 50
# RD_A2_A1is0
RD.cont_A2_A1is0 <- E.Y01 - E.Y00 # [1] 40
# RD_A2_A1is1
RD.cont_A2_A1is1 <- E.Y11 - E.Y10 # [1] 70

# interaction effects
# a.INT
a.INT.cont <- E.Y11 - E.Y10 - E.Y01 + E.Y00 # [1] 30


################################################################################
### 1) Simulations with binary outcomes
################################################################################
n.simu <- 100
library(MargIntTmle)
Q_form_Ybin <- c(Y.bin="Q.kplus1 ~ conf1 + conf2 + conf3 + fact.A1 * fact.A2")
g_form <- c("fact.A1 ~ conf1 + conf2",
            "fact.A2 ~ conf1 + conf3")

estim.gcomp.bin <- data.frame(k = 1:n.simu,
                              p00 = rep(NA, n.simu), p00.se = rep(NA, n.simu), p00.lb = rep(NA, n.simu), p00.ub = rep(NA, n.simu),
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.ub = rep(NA, n.simu),
                              p01 = rep(NA, n.simu), p01.se = rep(NA, n.simu), p01.lb = rep(NA, n.simu), p01.ub = rep(NA, n.simu),
                              p11 = rep(NA, n.simu), p11.se = rep(NA, n.simu), p11.lb = rep(NA, n.simu), p11.ub = rep(NA, n.simu),
                              RD_A1_A2is0=rep(NA, n.simu), RD_A1_A2is0.se=rep(NA, n.simu), RD_A1_A2is0.lb=rep(NA, n.simu), RD_A1_A2is0.ub=rep(NA, n.simu),
                              RD_A1_A2is1=rep(NA, n.simu), RD_A1_A2is1.se=rep(NA, n.simu), RD_A1_A2is1.lb=rep(NA, n.simu), RD_A1_A2is1.ub=rep(NA, n.simu),
                              RD_A2_A1is0=rep(NA, n.simu), RD_A2_A1is0.se=rep(NA, n.simu), RD_A2_A1is0.lb=rep(NA, n.simu), RD_A2_A1is0.ub=rep(NA, n.simu),
                              RD_A2_A1is1=rep(NA, n.simu), RD_A2_A1is1.se=rep(NA, n.simu), RD_A2_A1is1.lb=rep(NA, n.simu), RD_A2_A1is1.ub=rep(NA, n.simu),
                              RR_A1_A2is0=rep(NA, n.simu), RR_A1_A2is0.se=rep(NA, n.simu), RR_A1_A2is0.lb=rep(NA, n.simu), RR_A1_A2is0.ub=rep(NA, n.simu),
                              RR_A1_A2is1=rep(NA, n.simu), RR_A1_A2is1.se=rep(NA, n.simu), RR_A1_A2is1.lb=rep(NA, n.simu), RR_A1_A2is1.ub=rep(NA, n.simu),
                              RR_A2_A1is0=rep(NA, n.simu), RR_A2_A1is0.se=rep(NA, n.simu), RR_A2_A1is0.lb=rep(NA, n.simu), RR_A2_A1is0.ub=rep(NA, n.simu),
                              RR_A2_A1is1=rep(NA, n.simu), RR_A2_A1is1.se=rep(NA, n.simu), RR_A2_A1is1.lb=rep(NA, n.simu), RR_A2_A1is1.ub=rep(NA, n.simu),
                              a.INT=rep(NA, n.simu), a.INT.se=rep(NA, n.simu), a.INT.lb=rep(NA, n.simu), a.INT.ub=rep(NA, n.simu),
                              m.INT=rep(NA, n.simu), m.INT.se=rep(NA, n.simu), m.INT.lb=rep(NA, n.simu), m.INT.ub=rep(NA, n.simu),
                              RERI=rep(NA, n.simu), RERI.se=rep(NA, n.simu), RERI.lb=rep(NA, n.simu), RERI.ub=rep(NA, n.simu))

estim.iptw.bin <- data.frame(k = 1:n.simu,
                              p00 = rep(NA, n.simu), p00.se = rep(NA, n.simu), p00.lb = rep(NA, n.simu), p00.ub = rep(NA, n.simu),
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.ub = rep(NA, n.simu),
                              p01 = rep(NA, n.simu), p01.se = rep(NA, n.simu), p01.lb = rep(NA, n.simu), p01.ub = rep(NA, n.simu),
                              p11 = rep(NA, n.simu), p11.se = rep(NA, n.simu), p11.lb = rep(NA, n.simu), p11.ub = rep(NA, n.simu),
                              RD_A1_A2is0=rep(NA, n.simu), RD_A1_A2is0.se=rep(NA, n.simu), RD_A1_A2is0.lb=rep(NA, n.simu), RD_A1_A2is0.ub=rep(NA, n.simu),
                              RD_A1_A2is1=rep(NA, n.simu), RD_A1_A2is1.se=rep(NA, n.simu), RD_A1_A2is1.lb=rep(NA, n.simu), RD_A1_A2is1.ub=rep(NA, n.simu),
                              RD_A2_A1is0=rep(NA, n.simu), RD_A2_A1is0.se=rep(NA, n.simu), RD_A2_A1is0.lb=rep(NA, n.simu), RD_A2_A1is0.ub=rep(NA, n.simu),
                              RD_A2_A1is1=rep(NA, n.simu), RD_A2_A1is1.se=rep(NA, n.simu), RD_A2_A1is1.lb=rep(NA, n.simu), RD_A2_A1is1.ub=rep(NA, n.simu),
                              RR_A1_A2is0=rep(NA, n.simu), RR_A1_A2is0.se=rep(NA, n.simu), RR_A1_A2is0.lb=rep(NA, n.simu), RR_A1_A2is0.ub=rep(NA, n.simu),
                              RR_A1_A2is1=rep(NA, n.simu), RR_A1_A2is1.se=rep(NA, n.simu), RR_A1_A2is1.lb=rep(NA, n.simu), RR_A1_A2is1.ub=rep(NA, n.simu),
                              RR_A2_A1is0=rep(NA, n.simu), RR_A2_A1is0.se=rep(NA, n.simu), RR_A2_A1is0.lb=rep(NA, n.simu), RR_A2_A1is0.ub=rep(NA, n.simu),
                              RR_A2_A1is1=rep(NA, n.simu), RR_A2_A1is1.se=rep(NA, n.simu), RR_A2_A1is1.lb=rep(NA, n.simu), RR_A2_A1is1.ub=rep(NA, n.simu),
                              a.INT=rep(NA, n.simu), a.INT.se=rep(NA, n.simu), a.INT.lb=rep(NA, n.simu), a.INT.ub=rep(NA, n.simu),
                              m.INT=rep(NA, n.simu), m.INT.se=rep(NA, n.simu), m.INT.lb=rep(NA, n.simu), m.INT.ub=rep(NA, n.simu),
                              RERI=rep(NA, n.simu), RERI.se=rep(NA, n.simu), RERI.lb=rep(NA, n.simu), RERI.ub=rep(NA, n.simu))

estim.tmle.bin <- data.frame(k = 1:n.simu,
                              p00 = rep(NA, n.simu), p00.se = rep(NA, n.simu), p00.lb = rep(NA, n.simu), p00.ub = rep(NA, n.simu),
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.ub = rep(NA, n.simu),
                              p01 = rep(NA, n.simu), p01.se = rep(NA, n.simu), p01.lb = rep(NA, n.simu), p01.ub = rep(NA, n.simu),
                              p11 = rep(NA, n.simu), p11.se = rep(NA, n.simu), p11.lb = rep(NA, n.simu), p11.ub = rep(NA, n.simu),
                              RD_A1_A2is0=rep(NA, n.simu), RD_A1_A2is0.se=rep(NA, n.simu), RD_A1_A2is0.lb=rep(NA, n.simu), RD_A1_A2is0.ub=rep(NA, n.simu),
                              RD_A1_A2is1=rep(NA, n.simu), RD_A1_A2is1.se=rep(NA, n.simu), RD_A1_A2is1.lb=rep(NA, n.simu), RD_A1_A2is1.ub=rep(NA, n.simu),
                              RD_A2_A1is0=rep(NA, n.simu), RD_A2_A1is0.se=rep(NA, n.simu), RD_A2_A1is0.lb=rep(NA, n.simu), RD_A2_A1is0.ub=rep(NA, n.simu),
                              RD_A2_A1is1=rep(NA, n.simu), RD_A2_A1is1.se=rep(NA, n.simu), RD_A2_A1is1.lb=rep(NA, n.simu), RD_A2_A1is1.ub=rep(NA, n.simu),
                              RR_A1_A2is0=rep(NA, n.simu), RR_A1_A2is0.se=rep(NA, n.simu), RR_A1_A2is0.lb=rep(NA, n.simu), RR_A1_A2is0.ub=rep(NA, n.simu),
                              RR_A1_A2is1=rep(NA, n.simu), RR_A1_A2is1.se=rep(NA, n.simu), RR_A1_A2is1.lb=rep(NA, n.simu), RR_A1_A2is1.ub=rep(NA, n.simu),
                              RR_A2_A1is0=rep(NA, n.simu), RR_A2_A1is0.se=rep(NA, n.simu), RR_A2_A1is0.lb=rep(NA, n.simu), RR_A2_A1is0.ub=rep(NA, n.simu),
                              RR_A2_A1is1=rep(NA, n.simu), RR_A2_A1is1.se=rep(NA, n.simu), RR_A2_A1is1.lb=rep(NA, n.simu), RR_A2_A1is1.ub=rep(NA, n.simu),
                              a.INT=rep(NA, n.simu), a.INT.se=rep(NA, n.simu), a.INT.lb=rep(NA, n.simu), a.INT.ub=rep(NA, n.simu),
                              m.INT=rep(NA, n.simu), m.INT.se=rep(NA, n.simu), m.INT.lb=rep(NA, n.simu), m.INT.ub=rep(NA, n.simu),
                              RERI=rep(NA, n.simu), RERI.se=rep(NA, n.simu), RERI.lb=rep(NA, n.simu), RERI.ub=rep(NA, n.simu))

set.seed(54321)
for (i in 1:n.simu) {
  print(paste0("Simulation n° ",i))
  df.simu <- gen.sim.data(N = 2000, do.A1 = NULL, do.A2 = NULL)
  Ybin.ltmle <- int.ltmleMSM(data = subset(df.simu, select = -c(Y.cont)),
                             Qform = Q_form_Ybin,
                             gform = g_form,
                             Anodes = c("fact.A1", "fact.A2"),
                             Lnodes = c("conf1", "conf2", "conf3"),
                             Ynodes = c("Y.bin"),
                             SL.library = "glm",
                             gcomp = FALSE,
                             iptw.only = FALSE,
                             survivalOutcome = FALSE,
                             variance.method = "ic")

  ## TMLE estimation
  results.bin.tmle <- estim.int.effects(Ybin.ltmle, estimator = "tmle")$int.r
  # save results
  # p
  estim.tmle.bin$p00[estim.tmle.bin$k == i] <- results.bin.tmle$p[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p00.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.p[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p00.lb[estim.tmle.bin$k == i] <- results.bin.tmle$p.lo[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p00.ub[estim.tmle.bin$k == i] <- results.bin.tmle$p.up[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 0]

  estim.tmle.bin$p10[estim.tmle.bin$k == i] <- results.bin.tmle$p[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p10.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.p[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p10.lb[estim.tmle.bin$k == i] <- results.bin.tmle$p.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$p10.ub[estim.tmle.bin$k == i] <- results.bin.tmle$p.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]

  estim.tmle.bin$p01[estim.tmle.bin$k == i] <- results.bin.tmle$p[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p01.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.p[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p01.lb[estim.tmle.bin$k == i] <- results.bin.tmle$p.lo[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p01.ub[estim.tmle.bin$k == i] <- results.bin.tmle$p.up[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]

  estim.tmle.bin$p11[estim.tmle.bin$k == i] <- results.bin.tmle$p[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p11.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.p[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p11.lb[estim.tmle.bin$k == i] <- results.bin.tmle$p.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$p11.ub[estim.tmle.bin$k == i] <- results.bin.tmle$p.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  # RD
  estim.tmle.bin$RD_A1_A2is0[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RD_A1_A2is0.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.RD.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RD_A1_A2is0.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RD_A1_A2is0.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]

  estim.tmle.bin$RD_A1_A2is1[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A1_A2is1.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.RD.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A1_A2is1.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A1_A2is1.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A1.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]

  estim.tmle.bin$RD_A2_A1is0[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is0.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.RD.A2[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is0.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2.lo[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is0.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2.up[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]

  estim.tmle.bin$RD_A2_A1is1[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is1.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.RD.A2[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is1.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RD_A2_A1is1.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RD.A2.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  # RR
  estim.tmle.bin$RR_A1_A2is0[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RR_A1_A2is0.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.lnRR.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RR_A1_A2is0.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]
  estim.tmle.bin$RR_A1_A2is0.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 0]

  estim.tmle.bin$RR_A1_A2is1[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A1_A2is1.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.lnRR.A1[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A1_A2is1.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A1_A2is1.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A1.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]

  estim.tmle.bin$RR_A2_A1is0[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is0.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.lnRR.A2[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is0.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2.lo[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is0.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2.up[results.bin.tmle$A1 == 0 & results.bin.tmle$A2 == 1]

  estim.tmle.bin$RR_A2_A1is1[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is1.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.lnRR.A2[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is1.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RR_A2_A1is1.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RR.A2.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  # a.INT
  estim.tmle.bin$a.INT[estim.tmle.bin$k == i] <- results.bin.tmle$a.INT[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$a.INT.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.a.INT[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$a.INT.lb[estim.tmle.bin$k == i] <- results.bin.tmle$a.INT.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$a.INT.ub[estim.tmle.bin$k == i] <- results.bin.tmle$a.INT.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  # RERI
  estim.tmle.bin$RERI[estim.tmle.bin$k == i] <- results.bin.tmle$RERI[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RERI.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.lnRERI[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RERI.lb[estim.tmle.bin$k == i] <- results.bin.tmle$RERI.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$RERI.ub[estim.tmle.bin$k == i] <- results.bin.tmle$RERI.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  # m.INT
  estim.tmle.bin$m.INT[estim.tmle.bin$k == i] <- results.bin.tmle$m.INT[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$m.INT.se[estim.tmle.bin$k == i] <- results.bin.tmle$sd.ln.m.INT[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$m.INT.lb[estim.tmle.bin$k == i] <- results.bin.tmle$m.INT.lo[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]
  estim.tmle.bin$m.INT.ub[estim.tmle.bin$k == i] <- results.bin.tmle$m.INT.up[results.bin.tmle$A1 == 1 & results.bin.tmle$A2 == 1]

  ## IPTW estimation
  results.bin.iptw <- estim.int.effects(Ybin.ltmle, estimator = "iptw")$int.r
  # save results
  # p
  estim.iptw.bin$p00[estim.iptw.bin$k == i] <- results.bin.iptw$p[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p00.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.p[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p00.lb[estim.iptw.bin$k == i] <- results.bin.iptw$p.lo[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p00.ub[estim.iptw.bin$k == i] <- results.bin.iptw$p.up[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 0]

  estim.iptw.bin$p10[estim.iptw.bin$k == i] <- results.bin.iptw$p[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p10.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.p[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p10.lb[estim.iptw.bin$k == i] <- results.bin.iptw$p.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$p10.ub[estim.iptw.bin$k == i] <- results.bin.iptw$p.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]

  estim.iptw.bin$p01[estim.iptw.bin$k == i] <- results.bin.iptw$p[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p01.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.p[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p01.lb[estim.iptw.bin$k == i] <- results.bin.iptw$p.lo[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p01.ub[estim.iptw.bin$k == i] <- results.bin.iptw$p.up[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]

  estim.iptw.bin$p11[estim.iptw.bin$k == i] <- results.bin.iptw$p[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p11.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.p[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p11.lb[estim.iptw.bin$k == i] <- results.bin.iptw$p.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$p11.ub[estim.iptw.bin$k == i] <- results.bin.iptw$p.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  # RD
  estim.iptw.bin$RD_A1_A2is0[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RD_A1_A2is0.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.RD.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RD_A1_A2is0.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RD_A1_A2is0.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]

  estim.iptw.bin$RD_A1_A2is1[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A1_A2is1.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.RD.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A1_A2is1.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A1_A2is1.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A1.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]

  estim.iptw.bin$RD_A2_A1is0[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is0.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.RD.A2[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is0.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2.lo[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is0.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2.up[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]

  estim.iptw.bin$RD_A2_A1is1[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is1.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.RD.A2[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is1.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RD_A2_A1is1.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RD.A2.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  # RR
  estim.iptw.bin$RR_A1_A2is0[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RR_A1_A2is0.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.lnRR.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RR_A1_A2is0.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]
  estim.iptw.bin$RR_A1_A2is0.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 0]

  estim.iptw.bin$RR_A1_A2is1[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A1_A2is1.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.lnRR.A1[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A1_A2is1.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A1_A2is1.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A1.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]

  estim.iptw.bin$RR_A2_A1is0[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is0.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.lnRR.A2[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is0.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2.lo[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is0.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2.up[results.bin.iptw$A1 == 0 & results.bin.iptw$A2 == 1]

  estim.iptw.bin$RR_A2_A1is1[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is1.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.lnRR.A2[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is1.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RR_A2_A1is1.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RR.A2.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  # a.INT
  estim.iptw.bin$a.INT[estim.iptw.bin$k == i] <- results.bin.iptw$a.INT[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$a.INT.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.a.INT[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$a.INT.lb[estim.iptw.bin$k == i] <- results.bin.iptw$a.INT.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$a.INT.ub[estim.iptw.bin$k == i] <- results.bin.iptw$a.INT.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  # RERI
  estim.iptw.bin$RERI[estim.iptw.bin$k == i] <- results.bin.iptw$RERI[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RERI.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.lnRERI[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RERI.lb[estim.iptw.bin$k == i] <- results.bin.iptw$RERI.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$RERI.ub[estim.iptw.bin$k == i] <- results.bin.iptw$RERI.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  # m.INT
  estim.iptw.bin$m.INT[estim.iptw.bin$k == i] <- results.bin.iptw$m.INT[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$m.INT.se[estim.iptw.bin$k == i] <- results.bin.iptw$sd.ln.m.INT[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$m.INT.lb[estim.iptw.bin$k == i] <- results.bin.iptw$m.INT.lo[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]
  estim.iptw.bin$m.INT.ub[estim.iptw.bin$k == i] <- results.bin.iptw$m.INT.up[results.bin.iptw$A1 == 1 & results.bin.iptw$A2 == 1]

  ## g-computation estimation
  Ybin.gcomp <- int.ltmleMSM(data = subset(df.simu, select = -c(Y.cont)),
                             Qform = Q_form_Ybin,
                             gform = g_form,
                             Anodes = c("fact.A1", "fact.A2"),
                             Lnodes = c("conf1", "conf2", "conf3"),
                             Ynodes = c("Y.bin"),
                             SL.library = list(Q=list("glm"),g=list("glm")),
                             gcomp = TRUE,
                             iptw.only = FALSE,
                             survivalOutcome = FALSE,
                             variance.method = "ic",
                             B = 200)
  results.bin.gcomp <- estim.int.effects(Ybin.gcomp, estimator = "gcomp")$int.r
  # save results
  # p
  estim.gcomp.bin$p00[estim.gcomp.bin$k == i] <- results.bin.gcomp$p[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p00.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.p[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p00.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.lo[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p00.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.up[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 0]

  estim.gcomp.bin$p10[estim.gcomp.bin$k == i] <- results.bin.gcomp$p[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p10.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.p[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p10.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$p10.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]

  estim.gcomp.bin$p01[estim.gcomp.bin$k == i] <- results.bin.gcomp$p[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p01.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.p[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p01.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.lo[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p01.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.up[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]

  estim.gcomp.bin$p11[estim.gcomp.bin$k == i] <- results.bin.gcomp$p[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p11.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.p[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p11.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$p11.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$p.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  # RD
  estim.gcomp.bin$RD_A1_A2is0[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RD_A1_A2is0.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.RD.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RD_A1_A2is0.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RD_A1_A2is0.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]

  estim.gcomp.bin$RD_A1_A2is1[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A1_A2is1.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.RD.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A1_A2is1.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A1_A2is1.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A1.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]

  estim.gcomp.bin$RD_A2_A1is0[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is0.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.RD.A2[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is0.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2.lo[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is0.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2.up[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]

  estim.gcomp.bin$RD_A2_A1is1[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is1.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.RD.A2[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is1.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RD_A2_A1is1.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RD.A2.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  # RR
  estim.gcomp.bin$RR_A1_A2is0[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RR_A1_A2is0.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.lnRR.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RR_A1_A2is0.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]
  estim.gcomp.bin$RR_A1_A2is0.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 0]

  estim.gcomp.bin$RR_A1_A2is1[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A1_A2is1.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.lnRR.A1[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A1_A2is1.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A1_A2is1.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A1.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]

  estim.gcomp.bin$RR_A2_A1is0[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is0.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.lnRR.A2[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is0.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2.lo[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is0.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2.up[results.bin.gcomp$A1 == 0 & results.bin.gcomp$A2 == 1]

  estim.gcomp.bin$RR_A2_A1is1[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is1.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.lnRR.A2[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is1.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RR_A2_A1is1.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RR.A2.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  # a.INT
  estim.gcomp.bin$a.INT[estim.gcomp.bin$k == i] <- results.bin.gcomp$a.INT[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$a.INT.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.a.INT[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$a.INT.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$a.INT.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$a.INT.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$a.INT.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  # RERI
  estim.gcomp.bin$RERI[estim.gcomp.bin$k == i] <- results.bin.gcomp$RERI[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RERI.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.lnRERI[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RERI.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$RERI.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$RERI.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$RERI.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  # m.INT
  estim.gcomp.bin$m.INT[estim.gcomp.bin$k == i] <- results.bin.gcomp$m.INT[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$m.INT.se[estim.gcomp.bin$k == i] <- results.bin.gcomp$sd.ln.m.INT[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$m.INT.lb[estim.gcomp.bin$k == i] <- results.bin.gcomp$m.INT.lo[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
  estim.gcomp.bin$m.INT.ub[estim.gcomp.bin$k == i] <- results.bin.gcomp$m.INT.up[results.bin.gcomp$A1 == 1 & results.bin.gcomp$A2 == 1]
}

# save simulation results
saveRDS(estim.gcomp.bin, file = "./docs/estim_gcomp_bin")
saveRDS(estim.iptw.bin, file = "./docs/estim_iptw_bin")
saveRDS(estim.tmle.bin, file = "./docs/estim_tmle_bin")

estim_gcomp_bin <- readRDS(file = "./docs/estim_gcomp_bin")
estim_iptw_bin <- readRDS(file = "./docs/estim_iptw_bin")
estim_tmle_bin <- readRDS(file = "./docs/estim_tmle_bin")

perf.gcomp <- data.frame(bias = rep(NA, 15),
                         var = rep(NA, 15),
                         std.bias = rep(NA, 15),
                         mse = rep(NA, 15),
                         av.estim.se = rep(NA, 15),
                         cov = rep(NA, 15))
row.names(perf.gcomp) <- c("p00","p10","p01","p11",
                           "RD_A1_A2is0","RD_A1_A2is1","RD_A2_A1is0","RD_A2_A1is1",
                           "logRR_A1_A2is0","logRR_A1_A2is1","logRR_A2_A1is0","logRR_A2_A1is1",
                           "a.INT","log.m.INT","logRERI")
perf.iptw <- perf.tmle <- perf.gcomp

## performance summary - gcomp
# p00
perf.gcomp$bias[which(row.names(perf.gcomp) == "p00")] <- mean(estim_gcomp_bin$p00) - p00
perf.gcomp$var[which(row.names(perf.gcomp) == "p00")] <- mean((estim_gcomp_bin$p00 - mean(estim_gcomp_bin$p00))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "p00")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "p00")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "p00")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "p00")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "p00")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "p00")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "p00")] <- sqrt(mean(estim_gcomp_bin$p00.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "p00")] <- mean(as.numeric(p00 >= estim_gcomp_bin$p00.lb) & (p00 <= estim_gcomp_bin$p00.ub))
# p10
perf.gcomp$bias[which(row.names(perf.gcomp) == "p10")] <- mean(estim_gcomp_bin$p10) - p10
perf.gcomp$var[which(row.names(perf.gcomp) == "p10")] <- mean((estim_gcomp_bin$p10 - mean(estim_gcomp_bin$p10))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "p10")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "p10")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "p10")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "p10")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "p10")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "p10")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "p10")] <- sqrt(mean(estim_gcomp_bin$p10.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "p10")] <- mean(as.numeric(p10 >= estim_gcomp_bin$p10.lb) & (p10 <= estim_gcomp_bin$p10.ub))
# p01
perf.gcomp$bias[which(row.names(perf.gcomp) == "p01")] <- mean(estim_gcomp_bin$p01) - p01
perf.gcomp$var[which(row.names(perf.gcomp) == "p01")] <- mean((estim_gcomp_bin$p01 - mean(estim_gcomp_bin$p01))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "p01")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "p01")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "p01")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "p01")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "p01")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "p01")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "p01")] <- sqrt(mean(estim_gcomp_bin$p01.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "p01")] <- mean(as.numeric(p01 >= estim_gcomp_bin$p01.lb) & (p01 <= estim_gcomp_bin$p01.ub))
# p11
perf.gcomp$bias[which(row.names(perf.gcomp) == "p11")] <- mean(estim_gcomp_bin$p11) - p11
perf.gcomp$var[which(row.names(perf.gcomp) == "p11")] <- mean((estim_gcomp_bin$p11 - mean(estim_gcomp_bin$p11))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "p11")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "p11")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "p11")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "p11")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "p11")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "p11")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "p11")] <- sqrt(mean(estim_gcomp_bin$p11.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "p11")] <- mean(as.numeric(p11 >= estim_gcomp_bin$p11.lb) & (p11 <= estim_gcomp_bin$p11.ub))
# RD_A1_A2is0
perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- mean(estim_gcomp_bin$RD_A1_A2is0) - RD.bin_A1_A2is0
perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- mean((estim_gcomp_bin$RD_A1_A2is0 - mean(estim_gcomp_bin$RD_A1_A2is0))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is0")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is0")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is0")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is0")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- sqrt(mean(estim_gcomp_bin$RD_A1_A2is0.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "RD_A1_A2is0")] <- mean(as.numeric(RD.bin_A1_A2is0 >= estim_gcomp_bin$RD_A1_A2is0.lb) & (RD.bin_A1_A2is0 <= estim_gcomp_bin$RD_A1_A2is0.ub))
# RD_A1_A2is1
perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- mean(estim_gcomp_bin$RD_A1_A2is1) - RD.bin_A1_A2is1
perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- mean((estim_gcomp_bin$RD_A1_A2is1 - mean(estim_gcomp_bin$RD_A1_A2is1))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is1")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is1")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A1_A2is1")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A1_A2is1")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- sqrt(mean(estim_gcomp_bin$RD_A1_A2is1.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "RD_A1_A2is1")] <- mean(as.numeric(RD.bin_A1_A2is1 >= estim_gcomp_bin$RD_A1_A2is1.lb) & (RD.bin_A1_A2is1 <= estim_gcomp_bin$RD_A1_A2is1.ub))
# RD_A2_A1is0
perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- mean(estim_gcomp_bin$RD_A2_A1is0) - RD.bin_A2_A1is0
perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- mean((estim_gcomp_bin$RD_A2_A1is0 - mean(estim_gcomp_bin$RD_A2_A1is0))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is0")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is0")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is0")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is0")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- sqrt(mean(estim_gcomp_bin$RD_A2_A1is0.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "RD_A2_A1is0")] <- mean(as.numeric(RD.bin_A2_A1is0 >= estim_gcomp_bin$RD_A2_A1is0.lb) & (RD.bin_A2_A1is0 <= estim_gcomp_bin$RD_A2_A1is0.ub))
# RD_A2_A1is1
perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- mean(estim_gcomp_bin$RD_A2_A1is1) - RD.bin_A2_A1is1
perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- mean((estim_gcomp_bin$RD_A2_A1is1 - mean(estim_gcomp_bin$RD_A2_A1is1))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is1")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is1")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "RD_A2_A1is1")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "RD_A2_A1is1")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- sqrt(mean(estim_gcomp_bin$RR_A1_A2is0.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "RD_A2_A1is1")] <- mean(as.numeric(RD.bin_A2_A1is1 >= estim_gcomp_bin$RD_A2_A1is1.lb) & (RD.bin_A2_A1is1 <= estim_gcomp_bin$RD_A2_A1is1.ub))
# logRR_A1_A2is0
perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- mean(log(estim_gcomp_bin$RR_A1_A2is0)) - log(RR.bin_A1_A2is0)
perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- mean((log(estim_gcomp_bin$RR_A1_A2is0) - mean(log(estim_gcomp_bin$RR_A1_A2is0)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is0")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is0")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- sqrt(mean(estim_gcomp_bin$RD_A2_A1is1.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "logRR_A1_A2is0")] <- mean(as.numeric(RR.bin_A1_A2is0 > estim_gcomp_bin$RR_A1_A2is0.lb) & (RR.bin_A1_A2is0 < estim_gcomp_bin$RR_A1_A2is0.ub))
# logRR_A1_A2is1
perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- mean(log(estim_gcomp_bin$RR_A1_A2is1)) - log(RR.bin_A1_A2is1)
perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- mean((log(estim_gcomp_bin$RR_A1_A2is1) - mean(log(estim_gcomp_bin$RR_A1_A2is1)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is1")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A1_A2is1")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- sqrt(mean(estim_gcomp_bin$RR_A1_A2is1.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "logRR_A1_A2is1")] <- mean(as.numeric(RR.bin_A1_A2is1 > estim_gcomp_bin$RR_A1_A2is1.lb) & (RR.bin_A1_A2is1 < estim_gcomp_bin$RR_A1_A2is1.ub))
# logRR_A2_A1is0
perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- mean(log(estim_gcomp_bin$RR_A2_A1is0)) - log(RR.bin_A2_A1is0)
perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- mean((log(estim_gcomp_bin$RR_A2_A1is0) - mean(log(estim_gcomp_bin$RR_A2_A1is0)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is0")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is0")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- sqrt(mean(estim_gcomp_bin$RR_A2_A1is0.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "logRR_A2_A1is0")] <- mean(as.numeric(RR.bin_A2_A1is0 > estim_gcomp_bin$RR_A2_A1is0.lb) & (RR.bin_A2_A1is0 < estim_gcomp_bin$RR_A2_A1is0.ub))
# logRR_A2_A1is1
perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- mean(log(estim_gcomp_bin$RR_A2_A1is1)) - log(RR.bin_A2_A1is1)
perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- mean((log(estim_gcomp_bin$RR_A2_A1is1) - mean(log(estim_gcomp_bin$RR_A2_A1is1)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is1")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "logRR_A2_A1is1")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- sqrt(mean(estim_gcomp_bin$RR_A2_A1is1.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "logRR_A2_A1is1")] <- mean(as.numeric(RR.bin_A2_A1is1 > estim_gcomp_bin$RR_A2_A1is1.lb) & (RR.bin_A2_A1is1 < estim_gcomp_bin$RR_A2_A1is1.ub))
# a.INT
perf.gcomp$bias[which(row.names(perf.gcomp) == "a.INT")] <- mean(estim_gcomp_bin$a.INT) - a.INT.bin
perf.gcomp$var[which(row.names(perf.gcomp) == "a.INT")] <- mean((estim_gcomp_bin$a.INT - mean(estim_gcomp_bin$a.INT))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "a.INT")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "a.INT")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "a.INT")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "a.INT")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "a.INT")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "a.INT")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "a.INT")] <- sqrt(mean(estim_gcomp_bin$a.INT.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "a.INT")] <- mean(as.numeric(a.INT.bin >= estim_gcomp_bin$a.INT.lb) & (a.INT.bin <= estim_gcomp_bin$a.INT.ub))
# log.m.INT
perf.gcomp$bias[which(row.names(perf.gcomp) == "log.m.INT")] <- mean(log(estim_gcomp_bin$m.INT)) - log(m.INT.bin)
perf.gcomp$var[which(row.names(perf.gcomp) == "log.m.INT")] <- mean((log(estim_gcomp_bin$m.INT) - mean(log(estim_gcomp_bin$m.INT)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "log.m.INT")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "log.m.INT")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "log.m.INT")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "log.m.INT")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "log.m.INT")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "log.m.INT")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "log.m.INT")] <- sqrt(mean(estim_gcomp_bin$m.INT.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "log.m.INT")] <- mean(as.numeric(m.INT.bin > estim_gcomp_bin$m.INT.lb) & (m.INT.bin < estim_gcomp_bin$m.INT.ub))
# logRERI
perf.gcomp$bias[which(row.names(perf.gcomp) == "logRERI")] <- mean(log(estim_gcomp_bin$RERI)) - log(RERI.bin)
perf.gcomp$var[which(row.names(perf.gcomp) == "logRERI")] <- mean((log(estim_gcomp_bin$RERI) - mean(log(estim_gcomp_bin$RERI)))^2) * (nrow(estim_gcomp_bin)) / (nrow(estim_gcomp_bin) - 1)
perf.gcomp$std.bias[which(row.names(perf.gcomp) == "logRERI")] <- perf.gcomp$bias[which(row.names(perf.gcomp) == "logRERI")] / sqrt(perf.gcomp$var[which(row.names(perf.gcomp) == "logRERI")])
perf.gcomp$mse[which(row.names(perf.gcomp) == "logRERI")] <- perf.gcomp$var[which(row.names(perf.gcomp) == "logRERI")] + (perf.gcomp$bias[which(row.names(perf.gcomp) == "logRERI")])^2
perf.gcomp$av.estim.se[which(row.names(perf.gcomp) == "logRERI")] <- sqrt(mean(estim_gcomp_bin$RERI.se^2))
perf.gcomp$cov[which(row.names(perf.gcomp) == "logRERI")] <- mean(as.numeric(RERI.bin > estim_gcomp_bin$RERI.lb) & (RERI.bin < estim_gcomp_bin$RERI.ub))

## performance summary - IPTW
# p00
perf.iptw$bias[which(row.names(perf.iptw) == "p00")] <- mean(estim_iptw_bin$p00) - p00
perf.iptw$var[which(row.names(perf.iptw) == "p00")] <- mean((estim_iptw_bin$p00 - mean(estim_iptw_bin$p00))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "p00")] <- perf.iptw$bias[which(row.names(perf.iptw) == "p00")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "p00")])
perf.iptw$mse[which(row.names(perf.iptw) == "p00")] <- perf.iptw$var[which(row.names(perf.iptw) == "p00")] + (perf.iptw$bias[which(row.names(perf.iptw) == "p00")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "p00")] <- sqrt(mean(estim_iptw_bin$p00.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "p00")] <- mean(as.numeric(p00 >= estim_iptw_bin$p00.lb) & (p00 <= estim_iptw_bin$p00.ub))
# p10
perf.iptw$bias[which(row.names(perf.iptw) == "p10")] <- mean(estim_iptw_bin$p10) - p10
perf.iptw$var[which(row.names(perf.iptw) == "p10")] <- mean((estim_iptw_bin$p10 - mean(estim_iptw_bin$p10))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "p10")] <- perf.iptw$bias[which(row.names(perf.iptw) == "p10")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "p10")])
perf.iptw$mse[which(row.names(perf.iptw) == "p10")] <- perf.iptw$var[which(row.names(perf.iptw) == "p10")] + (perf.iptw$bias[which(row.names(perf.iptw) == "p10")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "p10")] <- sqrt(mean(estim_iptw_bin$p10.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "p10")] <- mean(as.numeric(p10 >= estim_iptw_bin$p10.lb) & (p10 <= estim_iptw_bin$p10.ub))
# p01
perf.iptw$bias[which(row.names(perf.iptw) == "p01")] <- mean(estim_iptw_bin$p01) - p01
perf.iptw$var[which(row.names(perf.iptw) == "p01")] <- mean((estim_iptw_bin$p01 - mean(estim_iptw_bin$p01))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "p01")] <- perf.iptw$bias[which(row.names(perf.iptw) == "p01")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "p01")])
perf.iptw$mse[which(row.names(perf.iptw) == "p01")] <- perf.iptw$var[which(row.names(perf.iptw) == "p01")] + (perf.iptw$bias[which(row.names(perf.iptw) == "p01")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "p01")] <- sqrt(mean(estim_iptw_bin$p01.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "p01")] <- mean(as.numeric(p01 >= estim_iptw_bin$p01.lb) & (p01 <= estim_iptw_bin$p01.ub))
# p11
perf.iptw$bias[which(row.names(perf.iptw) == "p11")] <- mean(estim_iptw_bin$p11) - p11
perf.iptw$var[which(row.names(perf.iptw) == "p11")] <- mean((estim_iptw_bin$p11 - mean(estim_iptw_bin$p11))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "p11")] <- perf.iptw$bias[which(row.names(perf.iptw) == "p11")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "p11")])
perf.iptw$mse[which(row.names(perf.iptw) == "p11")] <- perf.iptw$var[which(row.names(perf.iptw) == "p11")] + (perf.iptw$bias[which(row.names(perf.iptw) == "p11")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "p11")] <- sqrt(mean(estim_iptw_bin$p11.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "p11")] <- mean(as.numeric(p11 >= estim_iptw_bin$p11.lb) & (p11 <= estim_iptw_bin$p11.ub))
# RD_A1_A2is0
perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- mean(estim_iptw_bin$RD_A1_A2is0) - RD.bin_A1_A2is0
perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- mean((estim_iptw_bin$RD_A1_A2is0 - mean(estim_iptw_bin$RD_A1_A2is0))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is0")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is0")])
perf.iptw$mse[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is0")] + (perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is0")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- sqrt(mean(estim_iptw_bin$RD_A1_A2is0.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "RD_A1_A2is0")] <- mean(as.numeric(RD.bin_A1_A2is0 >= estim_iptw_bin$RD_A1_A2is0.lb) & (RD.bin_A1_A2is0 <= estim_iptw_bin$RD_A1_A2is0.ub))
# RD_A1_A2is1
perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- mean(estim_iptw_bin$RD_A1_A2is1) - RD.bin_A1_A2is1
perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- mean((estim_iptw_bin$RD_A1_A2is1 - mean(estim_iptw_bin$RD_A1_A2is1))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is1")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is1")])
perf.iptw$mse[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- perf.iptw$var[which(row.names(perf.iptw) == "RD_A1_A2is1")] + (perf.iptw$bias[which(row.names(perf.iptw) == "RD_A1_A2is1")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- sqrt(mean(estim_iptw_bin$RD_A1_A2is1.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "RD_A1_A2is1")] <- mean(as.numeric(RD.bin_A1_A2is1 >= estim_iptw_bin$RD_A1_A2is1.lb) & (RD.bin_A1_A2is1 <= estim_iptw_bin$RD_A1_A2is1.ub))
# RD_A2_A1is0
perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- mean(estim_iptw_bin$RD_A2_A1is0) - RD.bin_A2_A1is0
perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- mean((estim_iptw_bin$RD_A2_A1is0 - mean(estim_iptw_bin$RD_A2_A1is0))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is0")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is0")])
perf.iptw$mse[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is0")] + (perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is0")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- sqrt(mean(estim_iptw_bin$RD_A2_A1is0.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "RD_A2_A1is0")] <- mean(as.numeric(RD.bin_A2_A1is0 >= estim_iptw_bin$RD_A2_A1is0.lb) & (RD.bin_A2_A1is0 <= estim_iptw_bin$RD_A2_A1is0.ub))
# RD_A2_A1is1
perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- mean(estim_iptw_bin$RD_A2_A1is1) - RD.bin_A2_A1is1
perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- mean((estim_iptw_bin$RD_A2_A1is1 - mean(estim_iptw_bin$RD_A2_A1is1))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is1")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is1")])
perf.iptw$mse[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- perf.iptw$var[which(row.names(perf.iptw) == "RD_A2_A1is1")] + (perf.iptw$bias[which(row.names(perf.iptw) == "RD_A2_A1is1")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- sqrt(mean(estim_iptw_bin$RR_A1_A2is0.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "RD_A2_A1is1")] <- mean(as.numeric(RD.bin_A2_A1is1 >= estim_iptw_bin$RD_A2_A1is1.lb) & (RD.bin_A2_A1is1 <= estim_iptw_bin$RD_A2_A1is1.ub))
# logRR_A1_A2is0
perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- mean(log(estim_iptw_bin$RR_A1_A2is0)) - log(RR.bin_A1_A2is0)
perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- mean((log(estim_iptw_bin$RR_A1_A2is0) - mean(log(estim_iptw_bin$RR_A1_A2is0)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is0")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is0")])
perf.iptw$mse[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is0")] + (perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is0")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- sqrt(mean(estim_iptw_bin$RD_A2_A1is1.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "logRR_A1_A2is0")] <- mean(as.numeric(RR.bin_A1_A2is0 > estim_iptw_bin$RR_A1_A2is0.lb) & (RR.bin_A1_A2is0 < estim_iptw_bin$RR_A1_A2is0.ub))
# logRR_A1_A2is1
perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- mean(log(estim_iptw_bin$RR_A1_A2is1)) - log(RR.bin_A1_A2is1)
perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- mean((log(estim_iptw_bin$RR_A1_A2is1) - mean(log(estim_iptw_bin$RR_A1_A2is1)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is1")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is1")])
perf.iptw$mse[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- perf.iptw$var[which(row.names(perf.iptw) == "logRR_A1_A2is1")] + (perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A1_A2is1")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- sqrt(mean(estim_iptw_bin$RR_A1_A2is1.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "logRR_A1_A2is1")] <- mean(as.numeric(RR.bin_A1_A2is1 > estim_iptw_bin$RR_A1_A2is1.lb) & (RR.bin_A1_A2is1 < estim_iptw_bin$RR_A1_A2is1.ub))
# logRR_A2_A1is0
perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- mean(log(estim_iptw_bin$RR_A2_A1is0)) - log(RR.bin_A2_A1is0)
perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- mean((log(estim_iptw_bin$RR_A2_A1is0) - mean(log(estim_iptw_bin$RR_A2_A1is0)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is0")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is0")])
perf.iptw$mse[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is0")] + (perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is0")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- sqrt(mean(estim_iptw_bin$RR_A2_A1is0.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "logRR_A2_A1is0")] <- mean(as.numeric(RR.bin_A2_A1is0 > estim_iptw_bin$RR_A2_A1is0.lb) & (RR.bin_A2_A1is0 < estim_iptw_bin$RR_A2_A1is0.ub))
# logRR_A2_A1is1
perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- mean(log(estim_iptw_bin$RR_A2_A1is1)) - log(RR.bin_A2_A1is1)
perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- mean((log(estim_iptw_bin$RR_A2_A1is1) - mean(log(estim_iptw_bin$RR_A2_A1is1)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is1")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is1")])
perf.iptw$mse[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- perf.iptw$var[which(row.names(perf.iptw) == "logRR_A2_A1is1")] + (perf.iptw$bias[which(row.names(perf.iptw) == "logRR_A2_A1is1")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- sqrt(mean(estim_iptw_bin$RR_A2_A1is1.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "logRR_A2_A1is1")] <- mean(as.numeric(RR.bin_A2_A1is1 > estim_iptw_bin$RR_A2_A1is1.lb) & (RR.bin_A2_A1is1 < estim_iptw_bin$RR_A2_A1is1.ub))
# a.INT
perf.iptw$bias[which(row.names(perf.iptw) == "a.INT")] <- mean(estim_iptw_bin$a.INT) - a.INT.bin
perf.iptw$var[which(row.names(perf.iptw) == "a.INT")] <- mean((estim_iptw_bin$a.INT - mean(estim_iptw_bin$a.INT))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "a.INT")] <- perf.iptw$bias[which(row.names(perf.iptw) == "a.INT")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "a.INT")])
perf.iptw$mse[which(row.names(perf.iptw) == "a.INT")] <- perf.iptw$var[which(row.names(perf.iptw) == "a.INT")] + (perf.iptw$bias[which(row.names(perf.iptw) == "a.INT")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "a.INT")] <- sqrt(mean(estim_iptw_bin$a.INT.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "a.INT")] <- mean(as.numeric(a.INT.bin >= estim_iptw_bin$a.INT.lb) & (a.INT.bin <= estim_iptw_bin$a.INT.ub))
# log.m.INT
perf.iptw$bias[which(row.names(perf.iptw) == "log.m.INT")] <- mean(log(estim_iptw_bin$m.INT)) - log(m.INT.bin)
perf.iptw$var[which(row.names(perf.iptw) == "log.m.INT")] <- mean((log(estim_iptw_bin$m.INT) - mean(log(estim_iptw_bin$m.INT)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "log.m.INT")] <- perf.iptw$bias[which(row.names(perf.iptw) == "log.m.INT")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "log.m.INT")])
perf.iptw$mse[which(row.names(perf.iptw) == "log.m.INT")] <- perf.iptw$var[which(row.names(perf.iptw) == "log.m.INT")] + (perf.iptw$bias[which(row.names(perf.iptw) == "log.m.INT")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "log.m.INT")] <- sqrt(mean(estim_iptw_bin$m.INT.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "log.m.INT")] <- mean(as.numeric(m.INT.bin > estim_iptw_bin$m.INT.lb) & (m.INT.bin < estim_iptw_bin$m.INT.ub))
# logRERI
perf.iptw$bias[which(row.names(perf.iptw) == "logRERI")] <- mean(log(estim_iptw_bin$RERI)) - log(RERI.bin)
perf.iptw$var[which(row.names(perf.iptw) == "logRERI")] <- mean((log(estim_iptw_bin$RERI) - mean(log(estim_iptw_bin$RERI)))^2) * (nrow(estim_iptw_bin)) / (nrow(estim_iptw_bin) - 1)
perf.iptw$std.bias[which(row.names(perf.iptw) == "logRERI")] <- perf.iptw$bias[which(row.names(perf.iptw) == "logRERI")] / sqrt(perf.iptw$var[which(row.names(perf.iptw) == "logRERI")])
perf.iptw$mse[which(row.names(perf.iptw) == "logRERI")] <- perf.iptw$var[which(row.names(perf.iptw) == "logRERI")] + (perf.iptw$bias[which(row.names(perf.iptw) == "logRERI")])^2
perf.iptw$av.estim.se[which(row.names(perf.iptw) == "logRERI")] <- sqrt(mean(estim_iptw_bin$RERI.se^2))
perf.iptw$cov[which(row.names(perf.iptw) == "logRERI")] <- mean(as.numeric(RERI.bin > estim_iptw_bin$RERI.lb) & (RERI.bin < estim_iptw_bin$RERI.ub))

## performance summary - TMLE
# p00
perf.tmle$bias[which(row.names(perf.tmle) == "p00")] <- mean(estim_tmle_bin$p00) - p00
perf.tmle$var[which(row.names(perf.tmle) == "p00")] <- mean((estim_tmle_bin$p00 - mean(estim_tmle_bin$p00))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "p00")] <- perf.tmle$bias[which(row.names(perf.tmle) == "p00")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "p00")])
perf.tmle$mse[which(row.names(perf.tmle) == "p00")] <- perf.tmle$var[which(row.names(perf.tmle) == "p00")] + (perf.tmle$bias[which(row.names(perf.tmle) == "p00")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "p00")] <- sqrt(mean(estim_tmle_bin$p00.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "p00")] <- mean(as.numeric(p00 >= estim_tmle_bin$p00.lb) & (p00 <= estim_tmle_bin$p00.ub))
# p10
perf.tmle$bias[which(row.names(perf.tmle) == "p10")] <- mean(estim_tmle_bin$p10) - p10
perf.tmle$var[which(row.names(perf.tmle) == "p10")] <- mean((estim_tmle_bin$p10 - mean(estim_tmle_bin$p10))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "p10")] <- perf.tmle$bias[which(row.names(perf.tmle) == "p10")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "p10")])
perf.tmle$mse[which(row.names(perf.tmle) == "p10")] <- perf.tmle$var[which(row.names(perf.tmle) == "p10")] + (perf.tmle$bias[which(row.names(perf.tmle) == "p10")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "p10")] <- sqrt(mean(estim_tmle_bin$p10.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "p10")] <- mean(as.numeric(p10 >= estim_tmle_bin$p10.lb) & (p10 <= estim_tmle_bin$p10.ub))
# p01
perf.tmle$bias[which(row.names(perf.tmle) == "p01")] <- mean(estim_tmle_bin$p01) - p01
perf.tmle$var[which(row.names(perf.tmle) == "p01")] <- mean((estim_tmle_bin$p01 - mean(estim_tmle_bin$p01))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "p01")] <- perf.tmle$bias[which(row.names(perf.tmle) == "p01")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "p01")])
perf.tmle$mse[which(row.names(perf.tmle) == "p01")] <- perf.tmle$var[which(row.names(perf.tmle) == "p01")] + (perf.tmle$bias[which(row.names(perf.tmle) == "p01")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "p01")] <- sqrt(mean(estim_tmle_bin$p01.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "p01")] <- mean(as.numeric(p01 >= estim_tmle_bin$p01.lb) & (p01 <= estim_tmle_bin$p01.ub))
# p11
perf.tmle$bias[which(row.names(perf.tmle) == "p11")] <- mean(estim_tmle_bin$p11) - p11
perf.tmle$var[which(row.names(perf.tmle) == "p11")] <- mean((estim_tmle_bin$p11 - mean(estim_tmle_bin$p11))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "p11")] <- perf.tmle$bias[which(row.names(perf.tmle) == "p11")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "p11")])
perf.tmle$mse[which(row.names(perf.tmle) == "p11")] <- perf.tmle$var[which(row.names(perf.tmle) == "p11")] + (perf.tmle$bias[which(row.names(perf.tmle) == "p11")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "p11")] <- sqrt(mean(estim_tmle_bin$p11.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "p11")] <- mean(as.numeric(p11 >= estim_tmle_bin$p11.lb) & (p11 <= estim_tmle_bin$p11.ub))
# RD_A1_A2is0
perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- mean(estim_tmle_bin$RD_A1_A2is0) - RD.bin_A1_A2is0
perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- mean((estim_tmle_bin$RD_A1_A2is0 - mean(estim_tmle_bin$RD_A1_A2is0))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is0")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is0")])
perf.tmle$mse[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is0")] + (perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is0")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- sqrt(mean(estim_tmle_bin$RD_A1_A2is0.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "RD_A1_A2is0")] <- mean(as.numeric(RD.bin_A1_A2is0 >= estim_tmle_bin$RD_A1_A2is0.lb) & (RD.bin_A1_A2is0 <= estim_tmle_bin$RD_A1_A2is0.ub))
# RD_A1_A2is1
perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- mean(estim_tmle_bin$RD_A1_A2is1) - RD.bin_A1_A2is1
perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- mean((estim_tmle_bin$RD_A1_A2is1 - mean(estim_tmle_bin$RD_A1_A2is1))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is1")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is1")])
perf.tmle$mse[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- perf.tmle$var[which(row.names(perf.tmle) == "RD_A1_A2is1")] + (perf.tmle$bias[which(row.names(perf.tmle) == "RD_A1_A2is1")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- sqrt(mean(estim_tmle_bin$RD_A1_A2is1.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "RD_A1_A2is1")] <- mean(as.numeric(RD.bin_A1_A2is1 >= estim_tmle_bin$RD_A1_A2is1.lb) & (RD.bin_A1_A2is1 <= estim_tmle_bin$RD_A1_A2is1.ub))
# RD_A2_A1is0
perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- mean(estim_tmle_bin$RD_A2_A1is0) - RD.bin_A2_A1is0
perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- mean((estim_tmle_bin$RD_A2_A1is0 - mean(estim_tmle_bin$RD_A2_A1is0))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is0")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is0")])
perf.tmle$mse[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is0")] + (perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is0")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- sqrt(mean(estim_tmle_bin$RD_A2_A1is0.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "RD_A2_A1is0")] <- mean(as.numeric(RD.bin_A2_A1is0 >= estim_tmle_bin$RD_A2_A1is0.lb) & (RD.bin_A2_A1is0 <= estim_tmle_bin$RD_A2_A1is0.ub))
# RD_A2_A1is1
perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- mean(estim_tmle_bin$RD_A2_A1is1) - RD.bin_A2_A1is1
perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- mean((estim_tmle_bin$RD_A2_A1is1 - mean(estim_tmle_bin$RD_A2_A1is1))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is1")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is1")])
perf.tmle$mse[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- perf.tmle$var[which(row.names(perf.tmle) == "RD_A2_A1is1")] + (perf.tmle$bias[which(row.names(perf.tmle) == "RD_A2_A1is1")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- sqrt(mean(estim_tmle_bin$RR_A1_A2is0.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "RD_A2_A1is1")] <- mean(as.numeric(RD.bin_A2_A1is1 >= estim_tmle_bin$RD_A2_A1is1.lb) & (RD.bin_A2_A1is1 <= estim_tmle_bin$RD_A2_A1is1.ub))
# logRR_A1_A2is0
perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- mean(log(estim_tmle_bin$RR_A1_A2is0)) - log(RR.bin_A1_A2is0)
perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- mean((log(estim_tmle_bin$RR_A1_A2is0) - mean(log(estim_tmle_bin$RR_A1_A2is0)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is0")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is0")])
perf.tmle$mse[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is0")] + (perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is0")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- sqrt(mean(estim_tmle_bin$RD_A2_A1is1.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "logRR_A1_A2is0")] <- mean(as.numeric(RR.bin_A1_A2is0 > estim_tmle_bin$RR_A1_A2is0.lb) & (RR.bin_A1_A2is0 < estim_tmle_bin$RR_A1_A2is0.ub))
# logRR_A1_A2is1
perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- mean(log(estim_tmle_bin$RR_A1_A2is1)) - log(RR.bin_A1_A2is1)
perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- mean((log(estim_tmle_bin$RR_A1_A2is1) - mean(log(estim_tmle_bin$RR_A1_A2is1)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is1")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is1")])
perf.tmle$mse[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- perf.tmle$var[which(row.names(perf.tmle) == "logRR_A1_A2is1")] + (perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A1_A2is1")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- sqrt(mean(estim_tmle_bin$RR_A1_A2is1.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "logRR_A1_A2is1")] <- mean(as.numeric(RR.bin_A1_A2is1 > estim_tmle_bin$RR_A1_A2is1.lb) & (RR.bin_A1_A2is1 < estim_tmle_bin$RR_A1_A2is1.ub))
# logRR_A2_A1is0
perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- mean(log(estim_tmle_bin$RR_A2_A1is0)) - log(RR.bin_A2_A1is0)
perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- mean((log(estim_tmle_bin$RR_A2_A1is0) - mean(log(estim_tmle_bin$RR_A2_A1is0)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is0")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is0")])
perf.tmle$mse[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is0")] + (perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is0")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- sqrt(mean(estim_tmle_bin$RR_A2_A1is0.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "logRR_A2_A1is0")] <- mean(as.numeric(RR.bin_A2_A1is0 > estim_tmle_bin$RR_A2_A1is0.lb) & (RR.bin_A2_A1is0 < estim_tmle_bin$RR_A2_A1is0.ub))
# logRR_A2_A1is1
perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- mean(log(estim_tmle_bin$RR_A2_A1is1)) - log(RR.bin_A2_A1is1)
perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- mean((log(estim_tmle_bin$RR_A2_A1is1) - mean(log(estim_tmle_bin$RR_A2_A1is1)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is1")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is1")])
perf.tmle$mse[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- perf.tmle$var[which(row.names(perf.tmle) == "logRR_A2_A1is1")] + (perf.tmle$bias[which(row.names(perf.tmle) == "logRR_A2_A1is1")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- sqrt(mean(estim_tmle_bin$RR_A2_A1is1.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "logRR_A2_A1is1")] <- mean(as.numeric(RR.bin_A2_A1is1 > estim_tmle_bin$RR_A2_A1is1.lb) & (RR.bin_A2_A1is1 < estim_tmle_bin$RR_A2_A1is1.ub))
# a.INT
perf.tmle$bias[which(row.names(perf.tmle) == "a.INT")] <- mean(estim_tmle_bin$a.INT) - a.INT.bin
perf.tmle$var[which(row.names(perf.tmle) == "a.INT")] <- mean((estim_tmle_bin$a.INT - mean(estim_tmle_bin$a.INT))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "a.INT")] <- perf.tmle$bias[which(row.names(perf.tmle) == "a.INT")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "a.INT")])
perf.tmle$mse[which(row.names(perf.tmle) == "a.INT")] <- perf.tmle$var[which(row.names(perf.tmle) == "a.INT")] + (perf.tmle$bias[which(row.names(perf.tmle) == "a.INT")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "a.INT")] <- sqrt(mean(estim_tmle_bin$a.INT.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "a.INT")] <- mean(as.numeric(a.INT.bin >= estim_tmle_bin$a.INT.lb) & (a.INT.bin <= estim_tmle_bin$a.INT.ub))
# log.m.INT
perf.tmle$bias[which(row.names(perf.tmle) == "log.m.INT")] <- mean(log(estim_tmle_bin$m.INT)) - log(m.INT.bin)
perf.tmle$var[which(row.names(perf.tmle) == "log.m.INT")] <- mean((log(estim_tmle_bin$m.INT) - mean(log(estim_tmle_bin$m.INT)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "log.m.INT")] <- perf.tmle$bias[which(row.names(perf.tmle) == "log.m.INT")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "log.m.INT")])
perf.tmle$mse[which(row.names(perf.tmle) == "log.m.INT")] <- perf.tmle$var[which(row.names(perf.tmle) == "log.m.INT")] + (perf.tmle$bias[which(row.names(perf.tmle) == "log.m.INT")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "log.m.INT")] <- sqrt(mean(estim_tmle_bin$m.INT.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "log.m.INT")] <- mean(as.numeric(m.INT.bin > estim_tmle_bin$m.INT.lb) & (m.INT.bin < estim_tmle_bin$m.INT.ub))
# logRERI
perf.tmle$bias[which(row.names(perf.tmle) == "logRERI")] <- mean(log(estim_tmle_bin$RERI)) - log(RERI.bin)
perf.tmle$var[which(row.names(perf.tmle) == "logRERI")] <- mean((log(estim_tmle_bin$RERI) - mean(log(estim_tmle_bin$RERI)))^2) * (nrow(estim_tmle_bin)) / (nrow(estim_tmle_bin) - 1)
perf.tmle$std.bias[which(row.names(perf.tmle) == "logRERI")] <- perf.tmle$bias[which(row.names(perf.tmle) == "logRERI")] / sqrt(perf.tmle$var[which(row.names(perf.tmle) == "logRERI")])
perf.tmle$mse[which(row.names(perf.tmle) == "logRERI")] <- perf.tmle$var[which(row.names(perf.tmle) == "logRERI")] + (perf.tmle$bias[which(row.names(perf.tmle) == "logRERI")])^2
perf.tmle$av.estim.se[which(row.names(perf.tmle) == "logRERI")] <- sqrt(mean(estim_tmle_bin$RERI.se^2))
perf.tmle$cov[which(row.names(perf.tmle) == "logRERI")] <- mean(as.numeric(RERI.bin > estim_tmle_bin$RERI.lb) & (RERI.bin < estim_tmle_bin$RERI.ub))

# save simulation results
saveRDS(perf.gcomp, file = "./docs/perf_gcomp_bin")
saveRDS(perf.iptw, file = "./docs/perf_iptw_bin")
saveRDS(perf.tmle, file = "./docs/perf_tmle_bin")


################################################################################
### 2) simulations with quantitative outcomes
################################################################################
Q_form_Ycont <-c(Y.cont="Q.kplus1 ~ conf1 + conf2 + conf3 + fact.A1 * fact.A2")
g_form <- c("fact.A1 ~ conf1 + conf2",
            "fact.A2 ~ conf1 + conf3")
