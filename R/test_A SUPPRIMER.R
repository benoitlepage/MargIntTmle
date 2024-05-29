# ---------------------------------------------------------------------------- #
### A1 * A2 interaction effect ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
set.seed(12345)
df <- generate.data.multcat(N = 1000, b = param.causal.model.multcat())
head(df)
df <- data.frame(df[,c("conf1","conf2","conf3")],
                 behav.2 = ifelse(df$behav == 2, 1, 0),
                 behav.3 = ifelse(df$behav == 3, 1, 0),
                 env.2 = ifelse(df$env == 2, 1, 0),
                 env.3 = ifelse(df$env == 3, 1, 0),
                 hlth.outcome = df$hlth.outcome)
head(df)

Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
g_formulas = c("behav.2 ~ conf1 + conf2",
               "behav.3 ~ conf1 + conf2 + behav.2",
               "env.2 ~ conf1 + conf3",
               "env.3 ~ conf1 + conf3 + env.2")
# working.msm <- "Y ~ behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3"
#
# summary.measures.reg <- array(NA, dim = c(9,8,1))
# colnames(summary.measures.reg) <- c("behav.2","behav.3","env.2","env.3","behav.2:env.2","behav.2:env.3","env.2:behav.3","env.3:behav.3")
# summary.measures.reg[1,,1] <- c(0,0,0,0,0,0,0,0)
# summary.measures.reg[2,,1] <- c(1,0,0,0,0,0,0,0)
# summary.measures.reg[3,,1] <- c(0,1,0,0,0,0,0,0)
# summary.measures.reg[4,,1] <- c(0,0,1,0,0,0,0,0)
# summary.measures.reg[5,,1] <- c(0,0,0,1,0,0,0,0)
# summary.measures.reg[6,,1] <- c(1,0,1,0,1,0,0,0)
# summary.measures.reg[7,,1] <- c(1,0,0,1,0,1,0,0)
# summary.measures.reg[8,,1] <- c(0,1,1,0,0,0,1,0)
# summary.measures.reg[9,,1] <- c(0,1,0,1,0,0,0,1)

summary.measures.reg <- array(NA, dim = c(9,4,1))
colnames(summary.measures.reg) <- c("behav.2","behav.3","env.2","env.3")
summary.measures.reg[1,,1] <- c(0,0,0,0)
summary.measures.reg[2,,1] <- c(1,0,0,0)
summary.measures.reg[3,,1] <- c(0,1,0,0)
summary.measures.reg[4,,1] <- c(0,0,1,0)
summary.measures.reg[5,,1] <- c(0,0,0,1)
summary.measures.reg[6,,1] <- c(1,0,1,0)
summary.measures.reg[7,,1] <- c(1,0,0,1)
summary.measures.reg[8,,1] <- c(0,1,1,0)
summary.measures.reg[9,,1] <- c(0,1,0,1)

regimes.MSM <- array(NA, dim = c(nrow(df), 4, 9))
regimes.MSM[,,1] <- matrix(summary.measures.reg[1,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,2] <- matrix(summary.measures.reg[2,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,3] <- matrix(summary.measures.reg[3,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,4] <- matrix(summary.measures.reg[4,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,5] <- matrix(summary.measures.reg[5,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,6] <- matrix(summary.measures.reg[6,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,7] <- matrix(summary.measures.reg[7,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,8] <- matrix(summary.measures.reg[8,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)
regimes.MSM[,,9] <- matrix(summary.measures.reg[9,1:4,1], byrow = TRUE, nrow = nrow(df), ncol = 4)


ltmle_MSM <- ltmle::ltmleMSM(data = df,
                             Anodes = c("behav.2","behav.3","env.2","env.3"),
                             Cnodes = NULL,
                             Lnodes = NULL,
                             Ynodes = c("hlth.outcome"),
                             survivalOutcome = FALSE,
                             Qform = Q_formulas,
                             gform = g_formulas,
                             gbounds = c(0.01, 1),
                             Yrange = NULL,
                             deterministic.g.function = NULL,
                             SL.library = list(Q="SL.glm",
                                               g="SL.glm"),
                             SL.cvControl = list(),
                             regimes = regimes.MSM, # instead of abar
                             working.msm = working.msm,
                             summary.measures = summary.measures.reg,
                             final.Ynodes = NULL,
                             stratify = FALSE,
                             msm.weights = "empirical",
                             estimate.time = FALSE,
                             gcomp = FALSE,
                             iptw.only = FALSE,
                             deterministic.Q.function = NULL,
                             variance.method = "ic",
                             observation.weights = NULL,
                             id = NULL)
ltmle_MSM$msm
# Coefficients:
#      S1       S2       S3       S4       S5       S6       S7       S8       S9
# -2.4530   1.9935   0.6661   1.3743   2.8799  -0.4855  -0.4860  -0.3558   0.7099
#
# Degrees of Freedom: 9000 Total (i.e. Null);  8991 Residual
# Null Deviance:	    401.3
# Residual Deviance: 4.312 	AIC: NA

ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.4530092     1.9935141     0.6661055     1.3742659     2.8798853    -0.4855180    -0.4860343    -0.3558174     0.7099318

ltmle_MSM$beta.iptw
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.4382946     1.9640679     0.6288329     1.3675612     2.8631821    -0.4000502    -0.4649394    -0.3368394     0.7288878

ltmle_MSM_gcomp <- ltmle::ltmleMSM(data = df,
                                   Anodes = c("behav.2","behav.3","env.2","env.3"),
                                   Cnodes = NULL,
                                   Lnodes = NULL,
                                   Ynodes = c("hlth.outcome"),
                                   survivalOutcome = FALSE,
                                   Qform = Q_formulas,
                                   gform = g_formulas,
                                   gbounds = c(0.01, 1),
                                   Yrange = NULL,
                                   deterministic.g.function = NULL,
                                   SL.library = list(Q="SL.glm",
                                                     g="SL.glm"),
                                   SL.cvControl = list(),
                                   regimes = regimes.MSM, # instead of abar
                                   working.msm = working.msm,
                                   summary.measures = summary.measures.reg,
                                   final.Ynodes = NULL,
                                   stratify = FALSE,
                                   msm.weights = "empirical",
                                   estimate.time = FALSE,
                                   gcomp = TRUE,
                                   iptw.only = FALSE,
                                   deterministic.Q.function = NULL,
                                   variance.method = "ic",
                                   observation.weights = NULL,
                                   id = NULL)
ltmle_MSM_gcomp$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
# -2.24953551    1.64412480    0.31726384    1.19246185    2.72647970   -0.08805957   -0.29731438    0.02274258    0.81655811


# ---------------------------------------------------------------------------- #
### Now A2 is a effect modifier V in the set of baseline confounders ----
# ---------------------------------------------------------------------------- #
rm(list=ls())
set.seed(12345)
df <- generate.data.multcat(N = 1000, b = param.causal.model.multcat())

head(df)
df <- data.frame(df[,c("conf1","conf2","conf3")],
                 env.2 = ifelse(df$env == 2, 1, 0),
                 env.3 = ifelse(df$env == 3, 1, 0),
                 behav.2 = ifelse(df$behav == 2, 1, 0),
                 behav.3 = ifelse(df$behav == 3, 1, 0),
                 hlth.outcome = df$hlth.outcome)
head(df)

Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
g_formulas = c("behav.2 ~ conf1 + conf2",
               "behav.3 ~ conf1 + conf2 + behav.2")
working.msm <- "Y ~ behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3"

summary.measures.reg <- array(NA, dim = c(3,2,1))
colnames(summary.measures.reg) <- c("behav.2","behav.3")
summary.measures.reg[1,,1] <- c(0,0)
summary.measures.reg[2,,1] <- c(1,0)
summary.measures.reg[3,,1] <- c(0,1)

regimes.MSM <- array(NA, dim = c(nrow(df), 2, 3))
regimes.MSM[,,1] <- matrix(summary.measures.reg[1,,1], byrow = TRUE, nrow = nrow(df), ncol = 2)
regimes.MSM[,,2] <- matrix(summary.measures.reg[2,,1], byrow = TRUE, nrow = nrow(df), ncol = 2)
regimes.MSM[,,3] <- matrix(summary.measures.reg[3,,1], byrow = TRUE, nrow = nrow(df), ncol = 2)



ltmle_MSM <- ltmle::ltmleMSM(data = df,
                             Anodes = c("behav.2","behav.3"),
                             Cnodes = ,
                             Lnodes = NULL,
                             Ynodes = c("hlth.outcome"),
                             survivalOutcome = FALSE,
                             Qform = Q_formulas,
                             gform = g_formulas,
                             gbounds = c(0.01, 1),
                             Yrange = NULL,
                             deterministic.g.function = NULL,
                             SL.library = list(Q="SL.glm",
                                               g="SL.glm"),
                             SL.cvControl = list(),
                             regimes = regimes.MSM, # instead of abar
                             working.msm = working.msm,
                             summary.measures = summary.measures.reg,
                             final.Ynodes = NULL,
                             stratify = FALSE,
                             msm.weights = "empirical",
                             estimate.time = FALSE,
                             gcomp = FALSE,
                             iptw.only = FALSE,
                             deterministic.Q.function = NULL,
                             variance.method = "ic",
                             observation.weights = NULL,
                             id = NULL)
ltmle_MSM$msm
ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2686445     1.6019431     0.5048356     1.2089985     2.8154917    -0.2657478    -0.1442461    -0.2512775     0.6891194

ltmle_MSM$beta.iptw
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2548074     1.5685120     0.4691007     1.1964501     2.7993874    -0.1839724    -0.1181825    -0.2189340     0.7091611

ltmle_MSM_gcomp <- ltmle::ltmleMSM(data = df,
                                   Anodes = c("behav.2","behav.3"),
                                   Cnodes = ,
                                   Lnodes = NULL,
                                   Ynodes = c("hlth.outcome"),
                                   survivalOutcome = FALSE,
                                   Qform = Q_formulas,
                                   gform = g_formulas,
                                   gbounds = c(0.01, 1),
                                   Yrange = NULL,
                                   deterministic.g.function = NULL,
                                   SL.library = list(Q="SL.glm",
                                                     g="SL.glm"),
                                   SL.cvControl = list(),
                                   regimes = regimes.MSM, # instead of abar
                                   working.msm = working.msm,
                                   summary.measures = summary.measures.reg,
                                   final.Ynodes = NULL,
                                   stratify = FALSE,
                                   msm.weights = "empirical",
                                   estimate.time = FALSE,
                                   gcomp = TRUE,
                                   iptw.only = FALSE,
                                   deterministic.Q.function = NULL,
                                   variance.method = "ic",
                                   observation.weights = NULL,
                                   id = NULL)
ltmle_MSM_gcomp$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2457797     1.6437496     0.3079359     1.1916748     2.7257042    -0.0876299    -0.2964491     0.0236186     0.8184235
