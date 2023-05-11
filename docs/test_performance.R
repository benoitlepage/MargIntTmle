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
set.seed(12345)
p00 <- gen.sim.data(N = 1e6, do.A1 = 0, do.A2 = 0)
p10 <- gen.sim.data(N = 1e6, do.A1 = 1, do.A2 = 0)
p01 <- gen.sim.data(N = 1e6, do.A1 = 0, do.A2 = 1)
p11 <- gen.sim.data(N = 1e6, do.A1 = 1, do.A2 = 1)

### For Y binary
p00$mean.Y.bin # [1] 0.03441
p10$mean.Y.bin # [1] 0.059572
p01$mean.Y.bin # [1] 0.071873
p11$mean.Y.bin # [1] 0.212607

# Risk Difference
# RD_A1_A2is0
p10$mean.Y.bin - p00$mean.Y.bin # [1] 0.025162
# RD_A1_A2is1
p11$mean.Y.bin - p01$mean.Y.bin # [1] 0.140734
# RD_A2_A1is0
p01$mean.Y.bin - p00$mean.Y.bin # [1] 0.037463
# RD_A2_A1is1
p11$mean.Y.bin - p10$mean.Y.bin # [1] 0.153035

# Relative risk
# RR_A1_A2is0
p10$mean.Y.bin / p00$mean.Y.bin # [1] 1.731241
# RR_A1_A2is1
p11$mean.Y.bin / p01$mean.Y.bin # [1] 2.958093
# RR_A2_A1is0
p01$mean.Y.bin / p00$mean.Y.bin # [1] 2.088724
# RR_A2_A1is1
p11$mean.Y.bin / p10$mean.Y.bin # [1] 3.568908

# interaction effects
# a.INT
p11$mean.Y.bin - p10$mean.Y.bin - p01$mean.Y.bin + p00$mean.Y.bin # [1] 0.115572
# m.INT
(p11$mean.Y.bin * p00$mean.Y.bin) / (p10$mean.Y.bin * p01$mean.Y.bin)  # [1] 1.708655
# RERI
(p11$mean.Y.bin - p10$mean.Y.bin - p01$mean.Y.bin + p00$mean.Y.bin) / p00$mean.Y.bin # [1] 3.358675

### For Y continuous
p00$mean.Y.cont # [1] 132.0172
p10$mean.Y.cont # [1] 152.0218
p01$mean.Y.cont # [1] 172.0286
p11$mean.Y.cont # [1] 221.992

# Risk Difference
# RD_A1_A2is0
p10$mean.Y.cont - p00$mean.Y.cont # [1] 20.00453
# RD_A1_A2is1
p11$mean.Y.cont - p01$mean.Y.cont # [1] 49.96336
# RD_A2_A1is0
p01$mean.Y.cont - p00$mean.Y.cont # [1] 40.01139
# RD_A2_A1is1
p11$mean.Y.cont - p10$mean.Y.cont # [1] 69.97023

# interaction effects
# a.INT
p11$mean.Y.cont - p10$mean.Y.cont - p01$mean.Y.cont + p00$mean.Y.cont # [1] 29.95883


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
n.simu <- 10
library(MargIntTmle)
Q_form_Ybin <- c(Y.bin="Q.kplus1 ~ conf1 + conf2 + conf3 + fact.A1 * fact.A2")
g_form <- c("fact.A1 ~ conf1 + conf2",
            "fact.A2 ~ conf1 + conf3")

estim.gcomp.bin <- data.frame(k = 1:n.simu,
                              p00 = rep(NA, n.simu), p00.se = rep(NA, n.simu), p00.lb = rep(NA, n.simu), p00.ub = rep(NA, n.simu),
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.up = rep(NA, n.simu),
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
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.up = rep(NA, n.simu),
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
                              p10 = rep(NA, n.simu), p10.se = rep(NA, n.simu), p10.lb = rep(NA, n.simu), p10.up = rep(NA, n.simu),
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

set.seed(132435)
for (i in 1:n.simu) {
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

  est.bin.tmle <- estim.int.effects(Ybin.ltmle, estimator = "tmle")
  est.bin.iptw <- estim.int.effects(Ybin.ltmle, estimator = "iptw")


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
  est.bin.gcomp <- estim.int.effects(Ybin.gcomp, estimator = "gcomp")
}



### 2) simulations with quantitative outcomes
Q_form_Ycont <-c(Y.cont="Q.kplus1 ~ conf1 + conf2 + conf3 + fact.A1 * fact.A2")
g_form <- c("fact.A1 ~ conf1 + conf2",
            "fact.A2 ~ conf1 + conf3")
