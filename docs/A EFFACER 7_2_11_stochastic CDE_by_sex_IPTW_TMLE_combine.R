rm(list=ls())

### Upload the MSM beta coefficients and covariance matrices
betas_tmle <- vector("list", 10)
vars_tmle <- vector("list", 10)
betas_iptw <- vector("list", 10)
vars_iptw <- vector("list", 10)

for(i in 1:10) {
betas_tmle[[i]] <- readRDS(file = paste0("NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/betas_tmle",i,".rds"))[[i]]
vars_tmle[[i]] <- readRDS(file = paste0("NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/vars_tmle",i,".rds"))[[i]]
betas_iptw[[i]] <- readRDS(file = paste0("NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/betas_iptw",i,".rds"))[[i]]
vars_iptw[[i]] <- readRDS(file = paste0("NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/vars_iptw",i,".rds"))[[i]]
}

### pool the 10 estimations with MIcombine()
library(mitools)
MSM_CDEstoch_tmle_combined_bysex <- MIcombine(betas_tmle,vars_tmle)
MSM_CDEstoch_iptw_combined_bysex <- MIcombine(betas_iptw,vars_iptw)

saveRDS(MSM_CDEstoch_tmle_combined_bysex, file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/MSM_CDEstoch_tmle_combined_bysex.rds")
saveRDS(MSM_CDEstoch_iptw_combined_bysex, file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/MSM_CDEstoch_iptw_combined_bysex.rds")

summary(MSM_CDEstoch_tmle_combined_bysex)
# Multiple imputation results:
#   MIcombine.default(betas_tmle, vars_tmle)
#                   results           se        (lower        upper) missInfo
# b0          -3.716925e+00 0.2267048621 -4.1659252706 -3.2679243064     29 %
# ACE1         2.920475e-01 0.4045624898 -0.5095867185  1.0936816430     30 %
# ACE2         9.197738e-01 0.4202415188  0.0920960361  1.7474515021     20 %
# M           -2.202363e-03 0.0131934408 -0.0288753148  0.0244705883     50 %
# M2           1.308967e-04 0.0002145103 -0.0003078550  0.0005696484     59 %
# sex          4.090293e-01 0.3524663396 -0.2967938768  1.1148524517     42 %
# ACE1.M       1.620665e-03 0.0233817106 -0.0458110203  0.0490523501     53 %
# ACE1.M2     -8.039698e-05 0.0004069189 -0.0009220436  0.0007612497     65 %
# ACE2.M      -7.542756e-03 0.0265526521 -0.0628345424  0.0477490305     69 %
# ACE2.M2     -4.378774e-05 0.0005184490 -0.0011599578  0.0010723824     84 %
# ACE1.sex     4.211860e-01 0.5572844645 -0.6792650706  1.5216370877     24 %
# ACE2.sex    -9.337698e-01 0.9166426224 -2.8509828205  0.9834431721     71 %
# M.sex        1.630615e-03 0.0205661757 -0.0409296272  0.0441908579     66 %
# M2.sex      -3.755203e-05 0.0003079123 -0.0006770904  0.0006019863     68 %
# ACE1.M.sex  -6.465999e-03 0.0299521240 -0.0672462876  0.0543142904     53 %
# ACE1.M2.sex -4.465526e-05 0.0004681134 -0.0010050333  0.0009157228     60 %
# ACE2.M.sex   3.059048e-02 0.0481126496 -0.0730558712  0.1342368376     84 %
# ACE2.M2.sex -3.265056e-04 0.0007456914 -0.0019419830  0.0012889719     86 %

MSM_CDEstoch_tmle_combined_bysex$coefficients
#            b0          ACE1          ACE2             M            M2           sex        ACE1.M       ACE1.M2        ACE2.M       ACE2.M2 
# -3.716925e+00  2.920475e-01  9.197738e-01 -2.202363e-03  1.308967e-04  4.090293e-01  1.620665e-03 -8.039698e-05 -7.542756e-03 -4.378774e-05 
#     ACE1.sex      ACE2.sex         M.sex        M2.sex    ACE1.M.sex   ACE1.M2.sex    ACE2.M.sex   ACE2.M2.sex 
# 4.211860e-01 -9.337698e-01  1.630615e-03 -3.755203e-05 -6.465999e-03 -4.465526e-05  3.059048e-02 -3.265056e-04

MSM_CDEstoch_tmle_combined_bysex$variance
#                        b0          ACE1          ACE2             M            M2           sex        ACE1:M       ACE1:M2        ACE2:M
# b0           5.139509e-02 -4.871461e-02 -4.404772e-02 -2.388347e-03  2.842692e-05 -5.966467e-02  3.120219e-03 -4.383694e-05  2.040582e-03
# ACE1        -4.871461e-02  1.636708e-01  4.135395e-02  2.189803e-03 -2.727877e-05  6.952131e-02 -6.009865e-03  5.571322e-05  2.392405e-04
# ACE2        -4.404772e-02  4.135395e-02  1.766029e-01  2.147812e-03 -2.544056e-05  4.601597e-02 -2.250210e-03  2.666234e-05 -5.420517e-03
# M           -2.388347e-03  2.189803e-03  2.147812e-03  1.740669e-04 -2.639296e-06  3.075099e-03 -2.172764e-04  3.687203e-06 -1.420547e-04
# M2           2.842692e-05 -2.727877e-05 -2.544056e-05 -2.639296e-06  4.601465e-08 -3.953983e-05  3.106735e-06 -5.868903e-08  1.872368e-06
# sex         -5.966467e-02  6.952131e-02  4.601597e-02  3.075099e-03 -3.953983e-05  1.242325e-01 -4.670516e-03  6.245124e-05 -2.247923e-03
# ACE1:M       3.120219e-03 -6.009865e-03 -2.250210e-03 -2.172764e-04  3.106735e-06 -4.670516e-03  5.467044e-04 -8.682258e-06  2.663967e-04
# ACE1:M2     -4.383694e-05  5.571322e-05  2.666234e-05  3.687203e-06 -5.868903e-08  6.245124e-05 -8.682258e-06  1.655830e-07 -5.800285e-06
# ACE2:M       2.040582e-03  2.392405e-04 -5.420517e-03 -1.420547e-04  1.872368e-06 -2.247923e-03  2.663967e-04 -5.800285e-06  7.050433e-04
# ACE2:M2     -2.032605e-05 -2.709846e-05  4.164860e-05  2.180374e-06 -3.662294e-08  2.712552e-05 -4.289619e-06  1.110929e-07 -1.256201e-05
# ACE1:sex     5.256048e-02 -1.580214e-01 -4.338181e-02 -2.646284e-03  3.294173e-05 -1.244826e-01  7.159436e-03 -7.693624e-05  1.663024e-03
# ACE2:sex     8.131976e-02 -8.484044e-02 -2.209641e-01 -5.644972e-03  9.065014e-05 -1.387904e-01  7.371320e-03 -1.300049e-04  8.881441e-03
# M:sex        3.089609e-03 -3.445320e-03 -2.414554e-03 -2.227464e-04  3.335084e-06 -6.254699e-03  3.324990e-04 -5.157563e-06  1.886175e-04
# M2:sex      -3.911795e-05  4.203049e-05  3.045491e-05  3.317360e-06 -5.475213e-08  7.878998e-05 -4.783106e-06  8.089748e-08 -2.911097e-06
# ACE1:M:sex  -3.207551e-03  5.963021e-03  2.065543e-03  2.201866e-04 -3.059508e-06  6.842789e-03 -5.822737e-04  8.999854e-06 -3.551435e-04
# ACE1:M2:sex  4.112166e-05 -5.208839e-05 -2.891424e-05 -3.511198e-06  5.596659e-08 -8.119477e-05  8.673716e-06 -1.611194e-07  7.297560e-06
# ACE2:M:sex  -3.667056e-03  1.339123e-03  8.792734e-03  3.163951e-04 -5.023123e-06  6.248960e-03 -5.324020e-04  1.101899e-05 -9.607447e-04
# ACE2:M2:sex  3.555721e-05  9.581387e-06 -9.171535e-05 -3.899458e-06  6.626772e-08 -7.284621e-05  7.217701e-06 -1.593392e-07  1.602727e-05
# ACE2:M2      ACE1:sex      ACE2:sex         M:sex        M2:sex    ACE1:M:sex   ACE1:M2:sex    ACE2:M:sex   ACE2:M2:sex
# b0          -2.032605e-05  5.256048e-02  8.131976e-02  3.089609e-03 -3.911795e-05 -3.207551e-03  4.112166e-05 -3.667056e-03  3.555721e-05
# ACE1        -2.709846e-05 -1.580214e-01 -8.484044e-02 -3.445320e-03  4.203049e-05  5.963021e-03 -5.208839e-05  1.339123e-03  9.581387e-06
# ACE2         4.164860e-05 -4.338181e-02 -2.209641e-01 -2.414554e-03  3.045491e-05  2.065543e-03 -2.891424e-05  8.792734e-03 -9.171535e-05
# M            2.180374e-06 -2.646284e-03 -5.644972e-03 -2.227464e-04  3.317360e-06  2.201866e-04 -3.511198e-06  3.163951e-04 -3.899458e-06
# M2          -3.662294e-08  3.294173e-05  9.065014e-05  3.335084e-06 -5.475213e-08 -3.059508e-06  5.596659e-08 -5.023123e-06  6.626772e-08
# sex          2.712552e-05 -1.244826e-01 -1.387904e-01 -6.254699e-03  7.878998e-05  6.842789e-03 -8.119477e-05  6.248960e-03 -7.284621e-05
# ACE1:M      -4.289619e-06  7.159436e-03  7.371320e-03  3.324990e-04 -4.783106e-06 -5.822737e-04  8.673716e-06 -5.324020e-04  7.217701e-06
# ACE1:M2      1.110929e-07 -7.693624e-05 -1.300049e-04 -5.157563e-06  8.089748e-08  8.999854e-06 -1.611194e-07  1.101899e-05 -1.593392e-07
# ACE2:M      -1.256201e-05  1.663024e-03  8.881441e-03  1.886175e-04 -2.911097e-06 -3.551435e-04  7.297560e-06 -9.607447e-04  1.602727e-05
# ACE2:M2      2.687894e-07 -2.282919e-05 -9.163596e-05 -3.073909e-06  5.504374e-08  6.659503e-06 -1.503919e-07  1.662014e-05 -3.211431e-07
# ACE1:sex    -2.282919e-05  3.105660e-01  7.049685e-02  6.043424e-03 -7.283909e-05 -1.215549e-02  1.148935e-04 -3.332737e-03  5.069414e-05
# ACE2:sex    -9.163596e-05  7.049685e-02  8.402337e-01  9.035018e-03 -1.373832e-04 -7.827057e-03  1.622396e-04 -3.458294e-02  3.571667e-04
# M:sex       -3.073909e-06  6.043424e-03  9.035018e-03  4.229676e-04 -6.087197e-06 -4.617063e-04  6.567041e-06 -5.318729e-04  7.253653e-06
# M2:sex       5.504374e-08 -7.283909e-05 -1.373832e-04 -6.087197e-06  9.480996e-08  6.529359e-06 -1.033470e-07  8.565450e-06 -1.240938e-07
# ACE1:M:sex   6.659503e-06 -1.215549e-02 -7.827057e-03 -4.617063e-04  6.529359e-06  8.971297e-04 -1.260580e-05  6.853793e-04 -1.111598e-05
# ACE1:M2:sex -1.503919e-07  1.148935e-04  1.622396e-04  6.567041e-06 -1.033470e-07 -1.260580e-05  2.191301e-07 -1.457064e-05  2.377891e-07
# ACE2:M:sex   1.662014e-05 -3.332737e-03 -3.458294e-02 -5.318729e-04  8.565450e-06  6.853793e-04 -1.457064e-05  2.314827e-03 -3.298830e-05
# ACE2:M2:sex -3.211431e-07  5.069414e-05  3.571667e-04  7.253653e-06 -1.240938e-07 -1.111598e-05  2.377891e-07 -3.298830e-05  5.560557e-07

MSM_CDEstoch_tmle_combined_bysex$df
#        b0        ACE1        ACE2           M          M2         sex      ACE1.M     ACE1.M2      ACE2.M     ACE2.M2    ACE1.sex    ACE2.sex 
# 116.41976   111.42129   249.21157    39.61792    28.95603    56.92002    35.75299    23.06401    20.55714    13.45905   162.53680    19.19855 
#    M.sex      M2.sex  ACE1.M.sex ACE1.M2.sex  ACE2.M.sex ACE2.M2.sex 
# 22.84607    21.44058    35.41986    27.06702    13.37401    12.65252

100 * MSM_CDEstoch_tmle_combined_bysex$missinfo
#       b0        ACE1        ACE2           M          M2         sex      ACE1.M     ACE1.M2      ACE2.M     ACE2.M2    ACE1.sex    ACE2.sex 
# 29.01315    29.67201    19.64594    50.11848    58.52026    41.77443    52.74397    65.34749    69.03921    83.98851    24.45516    71.30882 
#    M.sex      M2.sex  ACE1.M.sex ACE1.M2.sex  ACE2.M.sex ACE2.M2.sex 
# 65.64601    67.67058    52.98940    60.47964    84.22789    86.34082


#_______________________________________________________________________________
#
# 4. Function to estimate probabilities and risk differences                ----
#_______________________________________________________________________________

Estimate_PSI <- function(msm) {
  coefs <- msm$coefficients
  var_IC <- msm$variance
  
  X.design <- function(beta,ace1,ace2,m,s) {
    X <- rep(NA, length(beta))
    names(X) <- names(beta)
    X["b0"] <- 1
    X["ACE1"] <- ace1
    X["ACE2"] <- ace2
    X["M"] <- m
    X["M2"] <- (m^2)
    X["sex"] <- s
    X["ACE1.M"] <- ace1 * m
    X["ACE1.M2"] <- ace1 * (m^2)
    X["ACE2.M"] <- ace2 * m
    X["ACE2.M2"] <- ace2 * (m^2)
    X["ACE1.sex"] <- ace1 * s 
    X["ACE2.sex"]  <- ace2 * s
    X["M.sex"] <- m * s
    X["M2.sex"] <- (m^2) * s
    X["ACE1.M.sex"] <- ace1 * m * s
    X["ACE1.M2.sex"] <- ace1 * (m^2) * s
    X["ACE2.M.sex"] <- ace2 * m * s
    X["ACE2.M2.sex"] <- ace2 * (m^2) * s
    return(X)
  }
  
  ## psi probabilities
  grid.p <- data.frame(ace1 = rep(c(0,1,0), 18), 
                       ace2 = rep(c(0,0,1),18), 
                       m = rep(seq(0,40, by = 5),each = 6), 
                       sex = rep(c(0,0,0,1,1,1),9))
  rownames(grid.p) <- paste0("p_ace", 1 * grid.p$ace1 + 2 * grid.p$ace2, "_m",grid.p$m,ifelse(grid.p$sex == 0,"w","m"))
  grid.p$p <- rep(NA, nrow(grid.p))
  for(r in 1:nrow(grid.p)) {
    grid.p$p[r] <- plogis(X.design(coefs, 
                                ace1 = grid.p$ace1[r], 
                                ace2 = grid.p$ace2[r], 
                                m = grid.p$m[r], 
                                s = grid.p$sex[r]) %*% as.matrix(coefs))
  }
  
  grid.RD <- data.frame(ace = rep(c(1,2),18), 
                        m = rep(seq(0,40, by = 5), each = 4), 
                        sex = rep(c(0,0,1,1), 9))
  rownames(grid.RD) <- paste0("RD_ace", grid.RD$ace, "v0_m",grid.RD$m,ifelse(grid.RD$sex == 0,"w","m"))
  grid.RD$RD <- rep(NA, nrow(grid.RD))
  for(r in 1:nrow(grid.RD)) {
    grid.RD$RD[r] <- (grid.p$p[grid.p$ace1 == ifelse(grid.RD$ace[r] == 1,1,0) & 
                                 grid.p$ace2 == ifelse(grid.RD$ace[r] == 2,1,0) & 
                                 grid.p$m == grid.RD$m[r] & 
                                 grid.p$sex == grid.RD$sex[r]] - 
                        grid.p$p[grid.p$ace1 == 0 & grid.p$ace2 == 0 & 
                                   grid.p$m == grid.RD$m[r] & 
                                   grid.p$sex == grid.RD$sex[r]])
  }

  ## gradient
  grad.p <- matrix(NA, nrow = nrow(grid.p), ncol = length(coefs))
  rownames(grad.p) <- rownames(grid.p)
  
  for(r in 1:nrow(grid.p)) {
    grad.p[r,] <- X.design(coefs, 
                           ace1 = grid.p$ace1[r], 
                           ace2 = grid.p$ace2[r], 
                           m = grid.p$m[r], 
                           s = grid.p$sex[r]) * grid.p$p[r] * (1 - grid.p$p[r])
  }
  
  grad.RD <- matrix(NA, nrow = nrow(grid.RD), ncol = length(coefs))
  rownames(grad.RD) <- rownames(grid.RD)
  for(r in 1:nrow(grid.RD)) {
    grad.RD[r,] <- (X.design(coefs, 
                           ace1 = ifelse(grid.RD$ace[r] == 1, 1,0), 
                           ace2 = ifelse(grid.RD$ace[r] == 2, 1,0), 
                           m = grid.RD$m[r], 
                           s = grid.RD$sex[r]) * 
                      grid.p$p[grid.p$ace1 == ifelse(grid.RD$ace[r] == 1,1,0) & 
                                 grid.p$ace2 == ifelse(grid.RD$ace[r] == 2,1,0) & 
                                 grid.p$m == grid.RD$m[r] & 
                                 grid.p$sex == grid.RD$sex[r]] * 
                      (1 - grid.p$p[grid.p$ace1 == ifelse(grid.RD$ace[r] == 1,1,0) & 
                                      grid.p$ace2 == ifelse(grid.RD$ace[r] == 2,1,0) & 
                                      grid.p$m == grid.RD$m[r] & 
                                      grid.p$sex == grid.RD$sex[r]])) - 
      (X.design(coefs, ace1 = 0, ace2 = 0, m = grid.RD$m[r], s = grid.RD$sex[r]) * 
         grid.p$p[grid.p$ace1 == 0 & grid.p$ace2 == 0 & 
                    grid.p$m == grid.RD$m[r] & 
                    grid.p$sex == grid.RD$sex[r]] * 
         (1 - grid.p$p[grid.p$ace1 == 0 & grid.p$ace2 == 0 & 
                         grid.p$m == grid.RD$m[r] & 
                         grid.p$sex == grid.RD$sex[r]]))
  }
  
  ## Standard errors
  grid.p$SE <- rep(NA, nrow(grid.p))
  for(r in 1:nrow(grid.p)) {
    grid.p$SE[r] <- sqrt(t(grad.p[r,]) %*% var_IC %*% grad.p[r,])
  }
  
  grid.RD$SE <- rep(NA, nrow(grid.RD))
  for(r in 1:nrow(grid.RD)) {
    grid.RD$SE[r] <- sqrt(t(grad.RD[r,]) %*% var_IC %*% grad.RD[r,])
  }
  
  ## Confidence intervals
  grid.p$low95CI <- grid.p$p - qnorm(0.975) * grid.p$SE
  grid.p$hi95CI <- grid.p$p + qnorm(0.975) * grid.p$SE
  
  grid.RD$low95CI <- grid.RD$RD - qnorm(0.975) * grid.RD$SE
  grid.RD$hi95CI <- grid.RD$RD + qnorm(0.975) * grid.RD$SE
  
  return(PSI = list(grid.p = grid.p,
                    grid.RD = grid.RD))
}

Psi_stoch_tmle <- Estimate_PSI(MSM_CDEstoch_tmle_combined_bysex)
Psi_stoch_iptw <- Estimate_PSI(MSM_CDEstoch_iptw_combined_bysex)
saveRDS(Psi_stoch_tmle, file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/Psi_stoch_tmle.rds")
saveRDS(Psi_stoch_iptw, file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/Psi_stoch_iptw.rds")



Psi_stoch_tmle$grid.p
#             ace1 ace2  m sex          p          SE      low95CI     hi95CI
# p_ace0_m0w     0    0  0   0 0.02373172 0.005252418  0.013437172 0.03402627
# p_ace1_m0w     1    0  0   0 0.03152697 0.010472277  0.011001682 0.05205225
# p_ace2_m0w     0    1  0   0 0.05747832 0.020263214  0.017763155 0.09719349
# p_ace0_m0m     0    0  0   1 0.03530132 0.008080349  0.019464125 0.05113851
# p_ace1_m0m     1    0  0   1 0.06948275 0.021662914  0.027024223 0.11194129
# p_ace2_m0m     0    1  0   1 0.03482777 0.024244679 -0.012690930 0.08234647  # negative lower bound => careful (ace = 2 in men)
# p_ace0_m5w     0    0  5   0 0.02355307 0.004155866  0.015407723 0.03169842
# p_ace1_m5w     1    0  5   0 0.03147675 0.009321352  0.013207235 0.04974626
# p_ace2_m5w     0    1  5   0 0.05500795 0.017792165  0.020135943 0.08987995
# p_ace0_m5m     0    0  5   1 0.03528344 0.006479171  0.022584497 0.04798238
# p_ace1_m5m     1    0  5   1 0.06770179 0.018336614  0.031762687 0.10364089
# p_ace2_m5m     0    1  5   1 0.03855189 0.022117933 -0.004798459 0.08190225 # negative lower bound => careful (ace = 2 in men)
# p_ace0_m10w    0    0 10   0 0.02352561 0.003399707  0.016862310 0.03018892
# p_ace1_m10w    1    0 10   0 0.03150356 0.008575003  0.014696861 0.04831025
# p_ace2_m10w    0    1 10   0 0.05285543 0.016751878  0.020022351 0.08568851
# p_ace0_m10m    0    0 10   1 0.03542470 0.005267314  0.025100956 0.04574845
# p_ace1_m10m    1    0 10   1 0.06586563 0.015651933  0.035188402 0.09654285
# p_ace2_m10m    0    1 10   1 0.04209472 0.019791216  0.003304647 0.08088479
# p_ace0_m15w    0    0 15   0 0.02364883 0.002940281  0.017885990 0.02941168
# p_ace1_m15w    1    0 15   0 0.03160758 0.008066938  0.015796674 0.04741849
# p_ace2_m15w    0    1 15   0 0.05099298 0.016273098  0.019098290 0.08288766
# p_ace0_m15m    0    0 15   1 0.03572695 0.004422883  0.027058259 0.04439564
# p_ace1_m15m    1    0 15   1 0.06398083 0.013488026  0.037544784 0.09041688
# p_ace2_m15m    0    1 15   1 0.04534434 0.017460722  0.011121955 0.07956673
# p_ace0_m20w    0    0 20   0 0.02392505 0.002699540  0.018634049 0.02921605
# p_ace1_m20w    1    0 20   0 0.03178956 0.007642198  0.016811132 0.04676800
# p_ace2_m20w    0    1 20   0 0.04939686 0.015771020  0.018486228 0.08030749
# p_ace0_m20m    0    0 20   1 0.03619414 0.003893160  0.028563689 0.04382460
# p_ace1_m20m    1    0 20   1 0.06205404 0.011727161  0.039069229 0.08503886
# p_ace2_m20m    0    1 20   1 0.04819287 0.015317421  0.018171272 0.07821446
# p_ace0_m25w    0    0 25   0 0.02435947 0.002581621  0.019299589 0.02941936
# p_ace1_m25w    1    0 25   0 0.03205080 0.007189886  0.017958884 0.04614272
# p_ace2_m25w    0    1 25   0 0.04804700 0.014968905  0.018708490 0.07738552
# p_ace0_m25m    0    0 25   1 0.03683243 0.003585380  0.029805213 0.04385964
# p_ace1_m25m    1    0 25   1 0.06009194 0.010282622  0.039938375 0.08024551
# p_ace2_m25m    0    1 25   1 0.05054208 0.013522373  0.024038717 0.07704544
# p_ace0_m30w    0    0 30   0 0.02496037 0.002514273  0.020032488 0.02988826
# p_ace1_m30w    1    0 30   0 0.03239316 0.006664707  0.019330577 0.04545575
# p_ace2_m30w    0    1 30   0 0.04692663 0.013826362  0.019827463 0.07402580
# p_ace0_m30m    0    0 30   1 0.03765027 0.003391714  0.031002637 0.04429791
# p_ace1_m30m    1    0 30   1 0.05810121 0.009130501  0.040205761 0.07599667
# p_ace2_m30m    0    1 30   1 0.05230871 0.012191890  0.028413045 0.07620437
# p_ace0_m35w    0    0 35   0 0.02573932 0.002484981  0.020868850 0.03060980
# p_ace1_m35w    1    0 35   0 0.03281911 0.006116904  0.020830198 0.04480802
# p_ace2_m35w    0    1 35   0 0.04602197 0.012526180  0.021471113 0.07057283
# p_ace0_m35m    0    0 35   1 0.03865863 0.003234197  0.032319726 0.04499754
# p_ace1_m35m    1    0 35   1 0.05608850 0.008336354  0.039749547 0.07242745
# p_ace2_m35m    0    1 35   1 0.05342890 0.011407573  0.031070468 0.07578733
# p_ace0_m40w    0    0 40   0 0.02671155 0.002569481  0.021675458 0.03174764
# p_ace1_m40w    1    0 40   0 0.03333173 0.005752671  0.022056701 0.04460676
# p_ace2_m40w    0    1 40   0 0.04532201 0.011563716  0.022657539 0.06798647
# p_ace0_m40m    0    0 40   1 0.03987117 0.003107352  0.033780875 0.04596147
# p_ace1_m40m    1    0 40   1 0.05406037 0.008051348  0.038280018 0.06984072
# p_ace2_m40m    0    1 40   1 0.05386165 0.011249261  0.031813502 0.07590979



Psi_stoch_tmle$grid.RD
#                ace  m sex            RD          SE       low95CI     hi95CI
# RD_ace1v0_m0w    1  0   0  0.0077952461 0.011552667 -0.0148475643 0.03043806
# RD_ace2v0_m0w    2  0   0  0.0337466017 0.020487598 -0.0064083527 0.07390156
# RD_ace1v0_m0m    1  0   1  0.0341814360 0.022621897 -0.0101566669 0.07851954
# RD_ace2v0_m0m    2  0   1 -0.0004735499 0.025520079 -0.0504919859 0.04954489
# RD_ace1v0_m5w    1  5   0  0.0079236786 0.009920738 -0.0115206108 0.02736797
# RD_ace2v0_m5w    2  5   0  0.0314548740 0.017918537 -0.0036648130 0.06657456
# RD_ace1v0_m5m    1  5   1  0.0324183496 0.019188074 -0.0051895843 0.07002628
# RD_ace2v0_m5m    2  5   1  0.0032684541 0.022894316 -0.0416035816 0.04814049
# RD_ace1v0_m10w   1 10   0  0.0079779438 0.008942301 -0.0095486449 0.02550453
# RD_ace2v0_m10w   2 10   0  0.0293298153 0.016778132 -0.0035547193 0.06221435
# RD_ace1v0_m10m   1 10   1  0.0304409235 0.016418719 -0.0017391749 0.06262102
# RD_ace2v0_m10m   2 10   1  0.0066700159 0.020229930 -0.0329799185 0.04631995
# RD_ace1v0_m15w   1 15   0  0.0079587474 0.008401192 -0.0085072862 0.02442478
# RD_ace2v0_m15w   2 15   0  0.0273441420 0.016225661 -0.0044575700 0.05914585
# RD_ace1v0_m15m   1 15   1  0.0282538789 0.014195805  0.0004306131 0.05607714
# RD_ace2v0_m15m   2 15   1  0.0096173909 0.017681729 -0.0250381611 0.04427294
# RD_ace1v0_m20w   1 20   0  0.0078645144 0.008064677 -0.0079419623 0.02367099
# RD_ace2v0_m20w   2 20   0  0.0254718104 0.015687774 -0.0052756618 0.05621928
# RD_ace1v0_m20m   1 20   1  0.0258598996 0.012403310  0.0015498591 0.05016994
# RD_ace2v0_m20m   2 20   1  0.0119987238 0.015407631 -0.0181996771 0.04219712
# RD_ace1v0_m25w   1 25   0  0.0076913281 0.007752336 -0.0075029708 0.02288563
# RD_ace2v0_m25w   2 25   0  0.0236875313 0.014894204 -0.0055045727 0.05287964
# RD_ace1v0_m25m   1 25   1  0.0232595156 0.010951628  0.0017947194 0.04472431
# RD_ace2v0_m25m   2 25   1  0.0137096525 0.013540647 -0.0128295282 0.04024883
# RD_ace1v0_m30w   1 30   0  0.0074327889 0.007375308 -0.0070225487 0.02188813
# RD_ace2v0_m30w   2 30   0  0.0219662601 0.013819018 -0.0051185174 0.04905104
# RD_ace1v0_m30m   1 30   1  0.0204509407 0.009805970  0.0012315925 0.03967029
# RD_ace2v0_m30m   2 30   1  0.0146584361 0.012177195 -0.0092084283 0.03852530
# RD_ace1v0_m35w   1 35   0  0.0070797872 0.006966640 -0.0065745771 0.02073415
# RD_ace2v0_m35w   2 35   0  0.0202826504 0.012676083 -0.0045620167 0.04512732
# RD_ace1v0_m35m   1 35   1  0.0174298654 0.009011351 -0.0002320585 0.03509179
# RD_ace2v0_m35m   2 35   1  0.0147702644 0.011399716 -0.0075727675 0.03711330
# RD_ace1v0_m40w   1 40   0  0.0066201808 0.006737950 -0.0065859589 0.01982632
# RD_ace2v0_m40w   2 40   0  0.0186104592 0.012000222 -0.0049095430 0.04213046
# RD_ace1v0_m40m   1 40   1  0.0141891979 0.008696458 -0.0028555456 0.03123394
# RD_ace2v0_m40m   2 40   1  0.0139904742 0.011324596 -0.0082053269 0.03618628

Psi_stoch_iptw <- Estimate_PSI(MSM_CDEstoch_iptw_combined_bysex)
Psi_stoch_iptw$grid.p
#             ace1 ace2  m sex          p          SE      low95CI     hi95CI
# p_ace0_m0w     0    0  0   0 0.02351156 0.005172009  0.013374607 0.03364851
# p_ace1_m0w     1    0  0   0 0.03067836 0.010714224  0.009678863 0.05167785
# p_ace2_m0w     0    1  0   0 0.05572453 0.020766565  0.015022812 0.09642625
# p_ace0_m0m     0    0  0   1 0.03841770 0.008599441  0.021563105 0.05527229
# p_ace1_m0m     1    0  0   1 0.07281419 0.023034179  0.027668032 0.11796036
# p_ace2_m0m     0    1  0   1 0.04447545 0.030334903 -0.014979865 0.10393077 # negative lower bound
# p_ace0_m5w     0    0  5   0 0.02432344 0.004166506  0.016157237 0.03248964
# p_ace1_m5w     1    0  5   0 0.03231710 0.009773132  0.013162110 0.05147208
# p_ace2_m5w     0    1  5   0 0.05570629 0.019004484  0.018458183 0.09295439
# p_ace0_m5m     0    0  5   1 0.03984666 0.007009665  0.026107968 0.05358535
# p_ace1_m5m     1    0  5   1 0.07334436 0.019640469  0.034849748 0.11183897
# p_ace2_m5m     0    1  5   1 0.04881603 0.027045501 -0.004192180 0.10182424 # negative lower bound
# p_ace0_m10w    0    0 10   0 0.02521629 0.003393810  0.018564546 0.03186804
# p_ace1_m10w    1    0 10   0 0.03390645 0.009091554  0.016087327 0.05172556
# p_ace2_m10w    0    1 10   0 0.05589036 0.018505923  0.019619415 0.09216130
# p_ace0_m10m    0    0 10   1 0.04128322 0.005690056  0.030130916 0.05243552
# p_ace1_m10m    1    0 10   1 0.07349449 0.016848602  0.040471836 0.10651714
# p_ace2_m10m    0    1 10   1 0.05301056 0.023709357  0.006541077 0.09948005
# p_ace0_m15w    0    0 15   0 0.02619673 0.002873091  0.020565580 0.03182789
# p_ace1_m15w    1    0 15   0 0.03543144 0.008591155  0.018593081 0.05226979
# p_ace2_m15w    0    1 15   0 0.05627863 0.018560422  0.019900874 0.09265639
# p_ace0_m15m    0    0 15   1 0.04272456 0.004697797  0.033517049 0.05193207
# p_ace1_m15m    1    0 15   1 0.07326243 0.014610216  0.044626931 0.10189793
# p_ace2_m15m    0    1 15   1 0.05695960 0.020568335  0.016646401 0.09727279
# p_ace0_m20w    0    0 20   0 0.02727215 0.002592846  0.022190263 0.03235403
# p_ace1_m20w    1    0 20   0 0.03687725 0.008154354  0.020895014 0.05285949
# p_ace2_m20w    0    1 20   0 0.05687511 0.018597397  0.020424876 0.09332533
# p_ace0_m20m    0    0 20   1 0.04416776 0.004065914  0.036198714 0.05213680
# p_ace1_m20m    1    0 20   1 0.07265151 0.012836018  0.047493375 0.09780964
# p_ace2_m20m    0    1 20   1 0.06056512 0.017825200  0.025628367 0.09550187
# p_ace0_m25w    0    0 25   0 0.02845077 0.002500787  0.023549314 0.03335222
# p_ace1_m25w    1    0 25   0 0.03822948 0.007664402  0.023207533 0.05325143
# p_ace2_m25w    0    1 25   0 0.05768593 0.018265693  0.021885827 0.09348603
# p_ace0_m25m    0    0 25   1 0.04560980 0.003763007  0.038234443 0.05298516
# p_ace1_m25m    1    0 25   1 0.07167048 0.011434131  0.049259995 0.09408096
# p_ace2_m25m    0    1 25   1 0.06373444 0.015598837  0.033161281 0.09430760
# p_ace0_m30w    0    0 30   0 0.02974176 0.002539991  0.024763473 0.03472005
# p_ace1_m30w    1    0 30   0 0.03947433 0.007050280  0.025656031 0.05329262
# p_ace2_m30w    0    1 30   0 0.05871951 0.017415473  0.024585810 0.09285321
# p_ace0_m30m    0    0 30   1 0.04704760 0.003689229  0.039816842 0.05427835
# p_ace1_m30m    1    0 30   1 0.07033334 0.010366597  0.050015181 0.09065150
# p_ace2_m30m    0    1 30   1 0.06638388 0.013908620  0.039123490 0.09364428
# p_ace0_m35w    0    0 35   0 0.03115536 0.002708983  0.025845856 0.03646487
# p_ace1_m35w    1    0 35   0 0.04059883 0.006344344  0.028164141 0.05303351
# p_ace2_m35w    0    1 35   0 0.05998664 0.016118580  0.028394799 0.09157847
# p_ace0_m35m    0    0 35   1 0.04847799 0.003742022  0.041143766 0.05581222
# p_ace1_m35m    1    0 35   1 0.06865906 0.009701291  0.049644879 0.08767324
# p_ace2_m35m    0    1 35   1 0.06844210 0.012737769  0.043476532 0.09340767
# p_ace0_m40w    0    0 40   0 0.03270295 0.003094860  0.026637139 0.03876877
# p_ace1_m40w    1    0 40   0 0.04159109 0.005785418  0.030251882 0.05293030
# p_ace2_m40w    0    1 40   0 0.06150062 0.014823740  0.032446627 0.09055462
# p_ace0_m40m    0    0 40   1 0.04989778 0.003892359  0.042268896 0.05752666
# p_ace1_m40m    1    0 40   1 0.06667121 0.009615681  0.047824824 0.08551760
# p_ace2_m40m    0    1 40   1 0.06985287 0.012183083  0.045974463 0.09373127

Psi_stoch_iptw$grid.RD
#                ace  m sex          RD          SE       low95CI     hi95CI
# RD_ace1v0_m0w    1  0   0 0.007166798 0.011768531 -0.0158990992 0.03023270
# RD_ace2v0_m0w    2  0   0 0.032212973 0.020951752 -0.0088517067 0.07327765
# RD_ace1v0_m0m    1  0   1 0.034396495 0.023989284 -0.0126216382 0.08141463
# RD_ace2v0_m0m    2  0   1 0.006057752 0.031927901 -0.0565197841 0.06863529
# RD_ace1v0_m5w    1  5   0 0.007993658 0.010412050 -0.0124135847 0.02840090
# RD_ace2v0_m5w    2  5   0 0.031382848 0.019129993 -0.0061112484 0.06887694
# RD_ace1v0_m5m    1  5   1 0.033497701 0.020539399 -0.0067587816 0.07375418
# RD_ace2v0_m5m    2  5   1 0.008969371 0.028007335 -0.0459239979 0.06386274
# RD_ace1v0_m10w   1 10   0 0.008690154 0.009496398 -0.0099224452 0.02730275
# RD_ace2v0_m10w   2 10   0 0.030674065 0.018519846 -0.0056241668 0.06697230
# RD_ace1v0_m10m   1 10   1 0.032211269 0.017680235 -0.0024413546 0.06686389
# RD_ace2v0_m10m   2 10   1 0.011727343 0.024145327 -0.0355966281 0.05905131
# RD_ace1v0_m15w   1 15   0 0.009234701 0.008939418 -0.0082862364 0.02675564
# RD_ace2v0_m15w   2 15   0 0.030081899 0.018462416 -0.0061037710 0.06626757
# RD_ace1v0_m15m   1 15   1 0.030537867 0.015361897  0.0004291028 0.06064663
# RD_ace2v0_m15m   2 15   1 0.014235034 0.020601918 -0.0261439826 0.05461405
# RD_ace1v0_m20w   1 20   0 0.009605107 0.008583964 -0.0072191533 0.02642937
# RD_ace2v0_m20w   2 20   0 0.029602958 0.018422103 -0.0065037005 0.06570962
# RD_ace1v0_m20m   1 20   1 0.028483750 0.013495500  0.0020330560 0.05493444
# RD_ace2v0_m20m   2 20   1 0.016397358 0.017608281 -0.0181142375 0.05090895
# RD_ace1v0_m25w   1 25   0 0.009778719 0.008258614 -0.0064078665 0.02596530
# RD_ace2v0_m25w   2 25   0 0.029235163 0.018068364 -0.0061781800 0.06464851
# RD_ace1v0_m25m   1 25   1 0.026060679 0.011990232  0.0025602561 0.04956110
# RD_ace2v0_m25m   2 25   1 0.018124639 0.015312029 -0.0118863875 0.04813566
# RD_ace1v0_m30w   1 30   0 0.009732563 0.007846220 -0.0056457458 0.02511087
# RD_ace2v0_m30w   2 30   0 0.028977747 0.017271057 -0.0048729026 0.06282840
# RD_ace1v0_m30m   1 30   1 0.023285741 0.010809441  0.0020996262 0.04447186
# RD_ace2v0_m30m   2 30   1 0.019336286 0.013743632 -0.0076007368 0.04627331
# RD_ace1v0_m35w   1 35   0 0.009443463 0.007350178 -0.0049626203 0.02384955
# RD_ace2v0_m35w   2 35   0 0.028831271 0.016131568 -0.0027860221 0.06044856
# RD_ace1v0_m35m   1 35   1 0.020181065 0.010028171  0.0005262108 0.03983592
# RD_ace2v0_m35m   2 35   1 0.019964106 0.012872100 -0.0052647459 0.04519296
# RD_ace1v0_m40w   1 40   0 0.008888140 0.006987647 -0.0048073956 0.02258368
# RD_ace2v0_m40w   2 40   0 0.028797671 0.015133790 -0.0008640134 0.05845935
# RD_ace1v0_m40m   1 40   1 0.016773431 0.009851394 -0.0025349462 0.03608181
# RD_ace2v0_m40m   2 40   1 0.019955087 0.012764818 -0.0050634961 0.04497367

#_______________________________________________________________________________
#
# 5. Figures                                                                ----
#_______________________________________________________________________________
library(ggplot2)

## 5.1 probabilities ----
#_______________________________________________________________________________
rm(list=ls())

MSM_CDEstoch_tmle_combined_bysex <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/MSM_CDEstoch_tmle_combined_bysex.rds")
MSM_CDEstoch_iptw_combined_bysex <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/MSM_CDEstoch_iptw_combined_bysex.rds")

Psi_stoch_tmle <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/Psi_stoch_tmle.rds")
Psi_stoch_iptw <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/stochastic_CDE/Psi_stoch_iptw.rds")

# import results of the ATE
PSI_ATE_iptw <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/results/PSI_ATE_iptw_bysex.rds")
PSI_ATE_tmle <- readRDS(file = "NCDS 58 CDE analysis 2022/working datasets/results/PSI_ATE_tmle_bysex.rds")

### 5.1.1 TMLE ----
#_______________________________________________________________________________

#### 5.1.1.a women ----
#_______________________________________________________________________________

# ACE = 0, women 
p.ace0.sex0 <- function(x) {
  ace1 = 0
  ace2 = 0
  sex = 0
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

# ACE = 1, women 
p.ace1.sex0 <- function(x) {
  ace1 = 1
  ace2 = 0
  sex = 0
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

# ACE = 2, women 
p.ace2.sex0 <- function(x) {
  ace1 = 0
  ace2 = 1
  sex = 0
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}


# do(ACE = 0)
g.ACE0.w <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace0.sex0, color = "red", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 0))$low95CI), 
            color = "red",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 0))$hi95CI), 
            color = "red",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE0w","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE0w","Estimate"]),
               linetype = "dotted", 
               color = "red") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE0w","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE0w","hi95CI"]), 
            linetype = "blank", 
            fill = "red",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.10)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=0), women")

# do(ACE = 1)
g.ACE1.w <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace1.sex0, color = "green4", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 0))$low95CI), 
            color = "green4",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 0))$hi95CI), 
            color = "green4",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE1w","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE1w","Estimate"]),
               linetype = "dotted", 
               color = "green4") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE1w","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE1w","hi95CI"]), 
            linetype = "blank", 
            fill = "green4",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.10)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=1), women")

# do(ACE = 2)
g.ACE2.w <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace2.sex0, color = "blue", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 0))$low95CI), 
            color = "blue",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 0))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 0))$hi95CI), 
            color = "blue",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE2w","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE2w","Estimate"]),
               linetype = "dotted", 
               color = "blue") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE2w","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE2w","hi95CI"]), 
            linetype = "blank", 
            fill = "blue",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.10)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=2), women")

#### 5.1.1.b men ----
#_______________________________________________________________________________

# ACE = 0, men 
p.ace0.sex1 <- function(x) {
  ace1 = 0
  ace2 = 0
  sex = 1
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

# ACE = 1, men 
p.ace1.sex1 <- function(x) {
  ace1 = 1
  ace2 = 0
  sex = 1
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

# ACE = 2, men 
p.ace2.sex1 <- function(x) {
  ace1 = 0
  ace2 = 1
  sex = 1
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           ace1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           ace2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           x * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           ace1 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           ace1 * (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           ace2 * x * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           ace2 * x^2 * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           ace1 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           ace2 * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           ace1 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           ace1 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           ace2 * x * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           ace2 * (x^2) * sex * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

# do(ACE = 0)
g.ACE0.m <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace0.sex1, color = "red", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 1))$low95CI), 
            color = "red",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 0 & sex == 1))$hi95CI), 
            color = "red",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE0m","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE0m","Estimate"]),
               linetype = "dotted", 
               color = "red") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE0m","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE0m","hi95CI"]), 
            linetype = "blank", 
            fill = "red",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.12)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=0), men")

# do(ACE = 1)
g.ACE1.m <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace1.sex1, color = "green4", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 1))$low95CI), 
            color = "green4",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 1 & ace2 == 0 & sex == 1))$hi95CI), 
            color = "green4",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE1m","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE1m","Estimate"]),
               linetype = "dotted", 
               color = "green4") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE1m","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE1m","hi95CI"]), 
            linetype = "blank", 
            fill = "green4",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.12)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=1), men")

# do(ACE = 2)
g.ACE2.m <- ggplot() + 
  # stochastic CDE
  geom_function(fun = p.ace2.sex1, color = "blue", linewidth = 1, linetype = "solid") +
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 1))$low95CI), 
            color = "blue",linetype = "dashed", linewidth = 1) + 
  geom_line(aes(x = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 1))$m, 
                y = subset(Psi_stoch_tmle$grid.p, subset = (ace1 == 0 & ace2 == 1 & sex == 1))$hi95CI), 
            color = "blue",linetype = "dashed", linewidth = 1) + 
  # ATE
  geom_segment(aes(x = 0, y = PSI_ATE_tmle["p_ACE2m","Estimate"],
                   xend = 40, yend = PSI_ATE_tmle["p_ACE2m","Estimate"]),
               linetype = "dotted", 
               color = "blue") + 
  geom_rect(aes(xmin = 0, xmax = 40, 
                ymin = PSI_ATE_tmle["p_ACE2m","low95CI"], 
                ymax = PSI_ATE_tmle["p_ACE2m","hi95CI"]), 
            linetype = "blank", 
            fill = "blue",
            alpha = 0.1) +
  xlim(c(0,40)) + ylim(c(0,0.12)) +
  xlab("P(M=1)") + ylab("P(Y=1)") + 
  ggtitle("Scenario do(ACE=2), men")

library(ggpubr)
ggarrange(g.ACE0.w, g.ACE1.w, g.ACE2.w, 
          g.ACE0.m, g.ACE1.m, g.ACE2.m,
          ncol = 3, nrow = 2)

## 5.2 risk differences ----
#_______________________________________________________________________________
### 5.2.1 Women ----
RD.ace1vs0.sex0 <- function(x) {
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           (1 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           (1 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           (1 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           (x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           ((x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           (1 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           (1 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"]) - 
    plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
             (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
             (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
             (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
             (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
             (x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
             ((x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
             (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
             (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
             (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
             (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

RD.ace2vs0.sex0 <- function(x) {
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           (1 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           (1 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           (1 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           (x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           ((x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           (1 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           (1 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"]) - 
    plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
             (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
             (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
             (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
             (0 * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
             (x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
             ((x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
             (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
             (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
             (0 * x * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
             (0 * (x^2) * 0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

ggplot() + 
  geom_function(fun = RD.ace1vs0.sex0, color = "green") +
  geom_function(fun = RD.ace2vs0.sex0, color = "blue") +
  xlim(c(0,40)) + 
  ggtitle("Risk difference according to P(M=1), in women, TMLE")

Psi_stoch_tmle$grid.RD$scenario <- rep("", nrow(Psi_stoch_tmle$grid.RD))
Psi_stoch_tmle$grid.RD$scenario[Psi_stoch_tmle$grid.RD$m %in% c(0,5,10,15,20,30,40)] <- "emulated to estimate MSM"
Psi_stoch_tmle$grid.RD$scenario[Psi_stoch_tmle$grid.RD$m == 0 & Psi_stoch_tmle$grid.RD$ace == 2] <- "extrapolated from MSM"

df1 <- subset(Psi_stoch_tmle$grid.RD, subset = (sex == 0))
df2 <- subset(Psi_stoch_tmle$grid.RD, subset = (sex == 0 & scenario != ""))
ggplot() + 
  geom_point(aes(x = df2$m, 
                 y = df2$RD, 
                 shape = as.factor(df2$scenario), 
                 color = as.factor(df2$ace)), size = 2) +
  geom_function(fun = RD.ace1vs0.sex0, color = "#F8766D") + # , linewidth = 0.8
  geom_function(fun = RD.ace2vs0.sex0, color = "#00BFC4") + # , linewidth = 0.8
  geom_line(aes(x = df1$m, y = df1$low95CI, color = as.factor(df1$ace)), linetype = "dashed") + 
  geom_line(aes(x = df1$m, y = df1$hi95CI, color = as.factor(df1$ace)), linetype = "dashed") + 
  xlab("Pr(M=1) in %") + ylab("Risk difference (in %)") +
  ggtitle("Risk difference according to P(M=1), in women, TMLE")


### 5.2.2 Men ----
RD.ace1vs0.sex1 <- function(x) {
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           (1 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           (1 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           (1 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           (x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           ((x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           (1 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           (1 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"]) - 
    plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
             (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
             (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
             (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
             (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
             (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
             (x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
             ((x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
             (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
             (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
             (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
             (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

RD.ace2vs0.sex1 <- function(x) {
  plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
           (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
           (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
           (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
           (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
           (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
           (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
           (1 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
           (1 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
           (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
           (1 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
           (x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
           ((x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
           (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
           (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
           (1 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
           (1 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"]) - 
    plogis(1 * MSM_CDEstoch_tmle_combined_bysex$coefficients["b0"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1"] + 
             (0) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2"] + 
             (x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M"] + 
             (x^2) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2"] + 
             (1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["sex"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2"] + 
             (0 * x) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M"] + 
             (0 * (x^2)) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2"] + 
             (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.sex"] + 
             (0 * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.sex"] + 
             (x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M.sex"] + 
             ((x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["M2.sex"] + 
             (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M.sex"] + 
             (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE1.M2.sex"] + 
             (0 * x * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M.sex"] + 
             (0 * (x^2) * 1) * MSM_CDEstoch_tmle_combined_bysex$coefficients["ACE2.M2.sex"])
}

ggplot() + 
  geom_function(fun = RD.ace1vs0.sex1, color = "green") +
  geom_function(fun = RD.ace2vs0.sex1, color = "blue") +
  xlim(c(0,40)) + 
  ggtitle("Risk difference according to P(M=1), in women, TMLE")