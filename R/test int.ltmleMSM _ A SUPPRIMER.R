

# ---------------------------------------------------------------------------- #
### interaction entre 2 expositions binaires ----
# ---------------------------------------------------------------------------- #
set.seed(12345)
df <- generate.data(N = 1000, b = param.causal.model())
head(df)

# specify the outcome model Q:
Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + sex * env")
# specify the treatment mechanisms g for each exposure A1 and A2:
g_formulas = c("sex ~ conf1 + conf2",
               "env ~ conf1 + conf3")

## choose a set of fitting libraries to pass to SuperLearner:
SL.library = list(Q=list("SL.glm"),g=list("SL.glm"))

## Apply the int.ltmleMSM function.
# In order to compute the TMLE and IPTW estimators, gcomp argument is set to FALSE.
interaction.ltmle <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("sex"),
                                  A2nodes = c("env"),
                                  Vnodes = NULL,
                                  Lnodes = c("conf1", "conf2", "conf3"),
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = FALSE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic")

interaction.ltmle$A1nodes
interaction.ltmle$A2nodes
interaction.ltmle$Vnodes
interaction.ltmle$Ynodes
interaction.ltmle$ltmle_MSM$beta
# (Intercept)         sex         env     sex:env
#  -2.3396212   1.7485790   0.8610451   1.8622763
interaction.ltmle$ltmle_MSM$msm
interaction.ltmle$ltmle_MSM$beta.iptw
# (Intercept)         sex         env     sex:env
#  -2.3225795   1.7459256   0.8358421   1.8627651
var(interaction.ltmle$ltmle_MSM$IC)
var(interaction.ltmle$ltmle_MSM$IC.iptw)

effects.tmle <- estim.int.effects(ltmle_MSM = interaction.ltmle,
                                  estimator = "tmle")

table_inter <- out.int.table(int.r = effects.tmle,
                             multipl = "RR",
                             labels.A1 = data.frame(effects.tmle$values.A1),
                             #   behav.2 behav.3 behav
                             # 1       0       0     1
                             # 2       1       0     2
                             # 3       0       1     3)
                             labels.A2 = data.frame(effects.tmle$values.A2),
                             #   env.2 env.3 env
                             # 1     0     0   1
                             # 2     1     0   2
                             # 3     0     1   3
                             labels.V = NULL)

library(kableExtra)
kbl(table_inter$out.table,
    caption = "Interaction effects estimated by TMLE") %>%
  kable_classic() %>%
  footnote(general = table_inter$interaction.effects)

# gcomp
interaction.gcomp <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("sex"),
                                  A2nodes = c("env"),
                                  Vnodes = NULL,
                                  Lnodes = c("conf1", "conf2", "conf3"),
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = TRUE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic",
                                  B = 100, # number of bootstrap samples
                                  boot.seed = 42) # seed for bootstrap
interaction.gcomp$A1nodes
interaction.gcomp$A2nodes
interaction.gcomp$Vnodes
interaction.gcomp$Ynodes
interaction.gcomp$ltmle_MSM$beta
# (Intercept)         sex         env     sex:env
#  -2.2766037   1.4812506   0.8825164   2.3446846
interaction.gcomp$ltmle_MSM$msm
boxplot(interaction.gcomp$bootstrap.res)
hist(interaction.gcomp$bootstrap.res$X.Intercept.)
hist(interaction.gcomp$bootstrap.res$sex)
hist(interaction.gcomp$bootstrap.res$env)
hist(interaction.gcomp$bootstrap.res$sex.env)
dim(interaction.gcomp$bootstrap.res)



# ---------------------------------------------------------------------------- #
### interaction entre 2 expositions à 3 classes ----
# ---------------------------------------------------------------------------- #
## Example 1
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

# Define Q and g formulas
# an A1 * A2 interaction term is recommended in the Q formula for the estimation
# of interaction effects
Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
g_formulas = c("behav.2 ~ conf1 + conf2",
               "behav.3 ~ conf1 + conf2 + behav.2",
               "env.2 ~ conf1 + conf3",
               "env.3 ~ conf1 + conf3 + env.2")

# Define SuperLearner libraries
SL.library = list(Q=list("SL.glm"),
                  g=list("SL.glm"))

# Estimate MSM parameters by IPTW and TMLE
interaction.ltmle <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("behav.2","behav.3"),
                                  A2nodes = c("env.2","env.3"),
                                  Vnodes = NULL,
                                  Lnodes = NULL,
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = FALSE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic")

interaction.ltmle$A1nodes
interaction.ltmle$A2nodes
interaction.ltmle$Vnodes
interaction.ltmle$Ynodes
interaction.ltmle$ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.4530092     1.9935141     0.6661055     1.3742659     2.8798853    -0.4855180    -0.4860343    -0.3558174     0.7099318
interaction.ltmle$ltmle_MSM$msm
interaction.ltmle$ltmle_MSM$beta.iptw
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.4382946     1.9640679     0.6288329     1.3675612     2.8631821    -0.4000502    -0.4649394    -0.3368394     0.7288878
var(interaction.ltmle$ltmle_MSM$IC)
var(interaction.ltmle$ltmle_MSM$IC.iptw)

effects.3cat.tmle <- estim.int.effects(ltmle_MSM = interaction.ltmle,
                                       estimator = "tmle")
effects.3cat.tmle$probs
effects.3cat.tmle$RD
effects.3cat.tmle$RR
effects.3cat.tmle$OR
effects.3cat.tmle$int




table_inter <- out.int.table(int.r = effects.3cat.tmle,
                             multipl = "RR",
                             labels.A1 = data.frame(effects.3cat.tmle$values.A1, behav = c(1,2,3)),
                             #   behav.2 behav.3 behav
                             # 1       0       0     1
                             # 2       1       0     2
                             # 3       0       1     3)
                             labels.A2 = data.frame(effects.3cat.tmle$values.A2, env = c(1,2,3)),
                             #   env.2 env.3 env
                             # 1     0     0   1
                             # 2     1     0   2
                             # 3     0     1   3
                             labels.V = NULL)

library(kableExtra)
kbl(table_inter$out.table,
    caption = "Interaction effects estimated by TMLE") %>%
  kable_classic() %>%
  footnote(general = table_inter$interaction.effects)

table_inter_OR <- out.int.table(int.r = effects.3cat.tmle,
                                multipl = "OR",
                                labels.A1 = data.frame(effects.3cat.tmle$values.A1, behav = c(1,2,3)),
                                labels.A2 = data.frame(effects.3cat.tmle$values.A2, env = c(1,2,3)),
                                labels.V = NULL)



effects.3cat.iptw <- estim.int.effects(ltmle_MSM = interaction.ltmle,
                                       estimator = "iptw")
effects.3cat.iptw$probs
effects.3cat.iptw$RD
effects.3cat.iptw$RR
effects.3cat.iptw$OR
effects.3cat.iptw$int

estim.int.effects(ltmle_MSM = interaction.ltmle,
                  estimator = "gcomp")
# ok cela indique une erreur



# Estimate MSM parameters by g-computation
interaction.gcomp <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("behav.2","behav.3"),
                                  A2nodes = c("env.2","env.3"),
                                  Vnodes = NULL,
                                  Lnodes = NULL,
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = TRUE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic",
                                  B = 100, # it should be at least 1000 or 2000
                                  boot.seed = 54321)
interaction.gcomp$A1nodes
interaction.gcomp$A2nodes
interaction.gcomp$Vnodes
interaction.gcomp$Ynodes
interaction.gcomp$ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
# -2.30462618    1.68915597    0.35220379    1.24696938    2.77394494   -0.12197048   -0.35760706   -0.01201771    0.79258807
interaction.gcomp$ltmle_MSM$msm
boxplot(interaction.gcomp$bootstrap.res)
hist(interaction.gcomp$bootstrap.res$X.Intercept.)
hist(interaction.gcomp$bootstrap.res$behav.2)
hist(interaction.gcomp$bootstrap.res$env.2)
hist(interaction.gcomp$bootstrap.res$env.3)
hist(interaction.gcomp$bootstrap.res$behav.3)
hist(interaction.gcomp$bootstrap.res$behav.2.env.2)
hist(interaction.gcomp$bootstrap.res$behav.2.env.3)
hist(interaction.gcomp$bootstrap.res$env.2.behav.3)
hist(interaction.gcomp$bootstrap.res$env.3.behav.3)
dim(interaction.gcomp$bootstrap.res)

effects.3cat.gcomp <- estim.int.effects(ltmle_MSM = interaction.gcomp,
                                       estimator = "gcomp")
effects.3cat.gcomp$probs
effects.3cat.gcomp$RD
effects.3cat.gcomp$RR
effects.3cat.gcomp$OR
effects.3cat.gcomp$int

## Example 2 - c(env.2, env3) are effect modifiers among baseline confounders
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

# Define Q and g formulas
# an (A1 * Effect.modifier) interaction term is recommended in the Q formula for the estimation
# of interaction effects
Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
g_formulas = c("behav.2 ~ conf1 + conf2",
               "behav.3 ~ conf1 + conf2 + behav.2") #c(behav2,behav.3) are the only exposure variables

# Define SuperLearner libraries
SL.library = list(Q = list("SL.glm"), g = list("SL.glm"))

# Estimate MSM parameters by IPTW and TMLE
interaction.ltmle <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("behav.2","behav.3"),
                                  A2nodes = NULL,
                                  Vnodes = c("env.2","env.3"),
                                  Lnodes = NULL,
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = FALSE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic")
interaction.ltmle$A1nodes
interaction.ltmle$A2nodes
interaction.ltmle$Vnodes
interaction.ltmle$Ynodes
interaction.ltmle$ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2686445     1.6019431     0.5048356     1.2089985     2.8154917    -0.2657478    -0.1442461    -0.2512775     0.6891194
interaction.ltmle$ltmle_MSM$msm
interaction.ltmle$ltmle_MSM$beta.iptw
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2548074     1.5685120     0.4691007     1.1964501     2.7993874    -0.1839724    -0.1181825    -0.2189340     0.7091611
var(interaction.ltmle$ltmle_MSM$IC)
var(interaction.ltmle$ltmle_MSM$IC.iptw)


interaction.gcomp <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  A1nodes = c("behav.2","behav.3"),
                                  A2nodes = NULL,
                                  Vnodes = c("env.2","env.3"),
                                  Lnodes = NULL,
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = TRUE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic",
                                  B = 100, # it should be at least 1000 or 2000
                                  boot.seed = 54321)
interaction.gcomp$A1nodes
interaction.gcomp$A2nodes
interaction.gcomp$Vnodes
interaction.gcomp$Ynodes
interaction.gcomp$ltmle_MSM$beta
# (Intercept)       behav.2         env.2         env.3       behav.3 behav.2:env.2 behav.2:env.3 env.2:behav.3 env.3:behav.3
#  -2.2457797     1.6437496     0.3079359     1.1916748     2.7257042    -0.0876299    -0.2964491     0.0236186     0.8184235
interaction.gcomp$ltmle_MSM$msm
boxplot(interaction.gcomp$bootstrap.res)
hist(interaction.gcomp$bootstrap.res$X.Intercept.)
hist(interaction.gcomp$bootstrap.res$behav.2)
hist(interaction.gcomp$bootstrap.res$env.2)
hist(interaction.gcomp$bootstrap.res$env.3)
hist(interaction.gcomp$bootstrap.res$behav.3)
hist(interaction.gcomp$bootstrap.res$behav.2.env.2)
hist(interaction.gcomp$bootstrap.res$behav.2.env.3)
hist(interaction.gcomp$bootstrap.res$env.2.behav.3)
hist(interaction.gcomp$bootstrap.res$env.3.behav.3)
dim(interaction.gcomp$bootstrap.res)

# OK seems to work !
