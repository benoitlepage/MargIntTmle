#' Estimate mean results, relative risks and differences
#'
#' @param ltmle_MSM an output from \code{int.ltmleMSM} function
#' @param estimator estimator of the marginal interaction effect. One of "gcomp", "iptw" or "tmle".
#'
#' @return \code{estim.int.effects} returns a list of 4 objects:
#'          \itemize{ \item \code{int.r} a data frame with the estimation of various quantities of interest for the interaction effect of A1 * A2 -> Y,
#'                    \item \code{Anodes} names of the exposure \code{c("A1","A2")},
#'                    \item \code{Ynodes} name of the outcome node,
#'                    \item \code{bootstrap.res} only for g-computation estimation: a data frame of length \code{int.ltmleMSM(B)} with the estimation of the MSM parameters from each bootstrap sample}
#' @export
#'
#' @examples
#' set.seed(12345)
#' df <- generate.data(N = 10000, b = param.causal.model())
#'
#' # Define Q and g formulas
#' # an A1 * A2 interaction term is recommended in the Q formula for the estimation
#' # of interaction effects
#' Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + sex * env")
#' g_formulas = c("sex ~ conf1 + conf2",
#'                "env ~ conf1 + conf3")
#'
#' # Define SuperLearner libraries
#' SL.library = list(Q=list("SL.glm", c("SL.glm", "screen.corP"),"SL.glmnet", "SL.mean"),
#'                   g=list("SL.glm", c("SL.glm", "screen.corP"),"SL.glmnet", "SL.mean"))
#' library(glmnet)
#'
#' # Estimate MSM parameters by IPTW and TMLE
#' interaction.ltmle <- int.ltmleMSM(data = df,
#'                                   Qform = Q_formulas,
#'                                   gform = g_formulas,
#'                                   Anodes = c("sex", "env"),
#'                                   Lnodes = c("conf1", "conf2", "conf3"),
#'                                   Ynodes = c("hlth.outcome"),
#'                                   SL.library = SL.library,
#'                                   gcomp = FALSE,
#'                                   iptw.only = FALSE,
#'                                   survivalOutcome = FALSE,
#'                                   variance.method = "ic")
#'
#' # Estimate quantities of interest for the interaction effect of (A1 * A2) on Y
#' estim.int.effects(interaction.ltmle, estimator = "tmle")
estim.int.effects <- function(ltmle_MSM = ltmle_MSM,
                              estimator = c("gcomp", "iptw", "tmle")) {

  data <- ltmle_MSM$data

  if(estimator == "gcomp") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == FALSE) stop("The ltmle function did not use the gcomp estimator, but the iptw +/- tmle estimator"))

    beta <- ltmle_MSM$ltmle_MSM$beta
  }

  if(estimator == "iptw") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == TRUE) stop("The ltmle function used the gcomp estimator, iptw is not available"))

    beta <- ltmle_MSM$ltmle_MSM$beta.iptw
    IC <- ltmle_MSM$ltmle_MSM$IC.iptw
  }

  if(estimator == "tmle") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == TRUE) stop("The ltmle function used the gcomp estimator, tmle is not available"))

    beta <- ltmle_MSM$ltmle_MSM$beta
    IC <- ltmle_MSM$ltmle_MSM$IC
  }

  # on va enregitrer l'ensemble des résultats pertinent dans une table de longueur k1 x k2
  int.r <- matrix(NA,
                  ncol = 34,
                  nrow = nlevels(as.factor(ltmle_MSM$data[,ltmle_MSM$Anodes[1]])) * nlevels(as.factor(ltmle_MSM$data[,ltmle_MSM$Anodes[2]])))
  int.r <- as.data.frame(int.r)
  names(int.r) <- c("A1","A2","p","sd.p","p.lo","p.up",
                    "RD.A1","sd.RD.A1","RD.A1.lo","RD.A1.up","RD.A2","sd.RD.A2","RD.A2.lo","RD.A2.up",
                    "RR.A1","sd.lnRR.A1","RR.A1.lo","RR.A1.up","RR.A2","sd.lnRR.A2","RR.A2.lo","RR.A2.up",
                    "a.INT", "sd.a.INT", "a.INT.lo", "a.INT.up","RERI","sd.RERI","RERI.lo","RERI.up",
                    "m.INT", "sd.ln.m.INT", "m.INT.lo", "m.INT.up" )
  int.r[,c("A1","A2")] <- expand.grid(c(0,1), c(0,1))

  # on peut retrouver les IC95% par delta method
  # A1 = 0 et A2 = 0
  int.r$p[int.r$A1 == 0 & int.r$A2 == 0] <- plogis(beta["(Intercept)"])

  # A1 = 1 et A2 = 0
  int.r$p[int.r$A1 == 1 & int.r$A2 == 0] <- plogis(beta["(Intercept)"] +
                                                     beta[ltmle_MSM$Anodes[1]]) # beta["A1"]

  # A1 = 0 et A2 = 1
  int.r$p[int.r$A1 == 0 & int.r$A2 == 1] <- plogis(beta["(Intercept)"] +
                                                     beta[ltmle_MSM$Anodes[2]]) # beta["A2"]

  # A1 = 1 et A2 = 1
  int.r$p[int.r$A1 == 1 & int.r$A2 == 1] <- plogis(beta["(Intercept)"] +
                                                     beta[ltmle_MSM$Anodes[1]] +         # beta["A1"]
                                                     beta[ltmle_MSM$Anodes[2]] +         # beta["A2"]
                                                     beta[paste0(ltmle_MSM$Anodes[1],":",ltmle_MSM$Anodes[2])])       # beta["A1:A2"]

  # RD.A1.A2is0
  int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 0] - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]

  # RD.A1.A2is1
  int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]

  # RD.A2.A1is0
  int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 1] - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]

  # RD.A2.A1is1
  int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]

  # RR.A1.A2is0
  int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- exp(log(int.r$p[int.r$A1 == 1 & int.r$A2 == 0]) - log(int.r$p[int.r$A1 == 0 & int.r$A2 == 0]))

  # RR.A1.A2is1
  int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) - log(int.r$p[int.r$A1 == 0 & int.r$A2 == 1]))

  # RR.A2.A1is0
  int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- exp(log(int.r$p[int.r$A1 == 0 & int.r$A2 == 1]) - log(int.r$p[int.r$A1 == 0 & int.r$A2 == 0]))

  # RR.A2.A1is1
  int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) - log(int.r$p[int.r$A1 == 1 & int.r$A2 == 0]))

  # additive interaction
  int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]

  # RERI
  int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] / int.r$p[int.r$A1 == 0 & int.r$A2 == 0]
    # exp(log(int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    #                                                      int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]) -
    #                                                  log(int.r$p[int.r$A1 == 0 & int.r$A2 == 0]))

  # multiplicative interaction
  int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) - log(int.r$p[int.r$A1 == 1 & int.r$A2 == 0]) -
                                                      log(int.r$p[int.r$A1 == 0 & int.r$A2 == 1]) + log(int.r$p[int.r$A1 == 0 & int.r$A2 == 0]))


  ## IC95%
  if(estimator == "iptw" | estimator == "tmle") {
    # A1 = 0 et A2 = 0
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),0,0,0)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0]
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0]

    # A1 = 1 et A2 = 0
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]),0,0)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0]
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0]

    # A1 = 0 et A2 = 1
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]), 0,
              int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]), 0)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1]
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1]

    # A1 = 1 et A2 = 1
    grad <- rep(int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]), 4)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1]

    # RD.A1.A2is0
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]), 0, 0)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0]
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0]

    # RD.A1.A2is1
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) )
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1]

    # RD.A2.A1is0
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]), 0,
              int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]), 0 )
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1]
    int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1]

    # RD.A2.A1is1
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]))
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1]

    # RR.A1.A2is0
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 0] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0], 0, 0)
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0])
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0])

    # RR.A1.A2is1
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              int.r$p[int.r$A1 == 0 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1] )
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                               qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                               qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1])

    # RR.A2.A1is0
    grad <- c(int.r$p[int.r$A1 == 0 & int.r$A2 == 0] - int.r$p[int.r$A1 == 0 & int.r$A2 == 1], 0,
              1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1], 0 )
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RR.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1])
    int.r$RR.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1])

    # RR.A2.A1is1
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 0] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              int.r$p[int.r$A1 == 1 & int.r$A2 == 0] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1])
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$RR.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$RR.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1])

    # additive interaction
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]) +
                int.r$p[int.r$A1 == 0 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
                int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]),
              int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) )
    grad.a.INT <- grad
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1]

    # RERI
    # sign.RERI <- ifelse(int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] < 0, -1, 1)
    # grad <- sign.RERI * c((int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
    #                          int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0]) -
    #                          int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1]) +
    #                          int.r$p[int.r$A1 == 0 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0])) /
    #                         (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    #                            int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]) -
    #                         (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),
    #                       (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
    #                          int.r$p[int.r$A1 == 1 & int.r$A2 == 0] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 0])) /
    #                         (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    #                            int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),
    #                       (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1]) -
    #                          int.r$p[int.r$A1 == 0 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 1])) /
    #                         (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    #                            int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]),
    #                       (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] * (1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1])) /
    #                         (int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
    #                            int.r$p[int.r$A1 == 0 & int.r$A2 == 1] + int.r$p[int.r$A1 == 0 & int.r$A2 == 0]) )
    # v <- t(grad) %*% var(IC) %*% grad
    # int.r$sd.lnRERI[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))
    #
    # abs.RERI <- sign.RERI * int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1]
    # abs.RERI.lo <- exp(log(abs.RERI) -
    #                       qnorm(0.975) * int.r$sd.lnRERI[int.r$A1 == 1 & int.r$A2 == 1])
    # abs.RERI.up <- exp(log(abs.RERI) +
    #                       qnorm(0.975) * int.r$sd.lnRERI[int.r$A1 == 1 & int.r$A2 == 1])
    #
    # int.r$RERI.lo[int.r$A1 == 1 & int.r$A2 == 1] <- sign.RERI * abs.RERI.lo
    # int.r$RERI.up[int.r$A1 == 1 & int.r$A2 == 1] <- sign.RERI * abs.RERI.up
    grad <- (grad.a.INT / int.r$p[int.r$A1 == 0 & int.r$A2 == 0]) -
      (int.r[int.r$A1 == 1 & int.r$A2 == 1,"a.INT"]) * (1 - int.r$p[int.r$A1 == 0 & int.r$A2 == 0]) / int.r$p[int.r$A1 == 0 & int.r$A2 == 0]
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data)) # sd
    int.r$RERI.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] - qnorm(0.975) * int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$RERI.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] + qnorm(0.975) * int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1]

    # multiplicative interaction
    grad <- c(int.r$p[int.r$A1 == 1 & int.r$A2 == 0] + int.r$p[int.r$A1 == 0 & int.r$A2 == 1] -
                int.r$p[int.r$A1 == 1 & int.r$A2 == 1] - int.r$p[int.r$A1 == 0 & int.r$A2 == 0],
              int.r$p[int.r$A1 == 1 & int.r$A2 == 0] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              int.r$p[int.r$A1 == 0 & int.r$A2 == 1] - int.r$p[int.r$A1 == 1 & int.r$A2 == 1],
              1 - int.r$p[int.r$A1 == 1 & int.r$A2 == 1])
    v <- t(grad) %*% var(IC) %*% grad
    int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- sqrt(v / nrow(ltmle_MSM$data))

    int.r$m.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$m.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1])

    bootstrap.res <- ltmle_MSM$bootstrap.res
  }

  if(estimator == "gcomp") {
    ltmle_MSM$bootstrap.res$p.A1_0.A2_0 <- plogis(ltmle_MSM$bootstrap.res$beta.Intercept)
    ltmle_MSM$bootstrap.res$p.A1_1.A2_0 <- plogis(ltmle_MSM$bootstrap.res$beta.Intercept +
                                                    ltmle_MSM$bootstrap.res$beta.A1)
    ltmle_MSM$bootstrap.res$p.A1_0.A2_1 <- plogis(ltmle_MSM$bootstrap.res$beta.Intercept +
                                                    ltmle_MSM$bootstrap.res$beta.A2)
    ltmle_MSM$bootstrap.res$p.A1_1.A2_1 <- plogis(ltmle_MSM$bootstrap.res$beta.Intercept +
                                                    ltmle_MSM$bootstrap.res$beta.A1 +
                                                    ltmle_MSM$bootstrap.res$beta.A2 +
                                                    ltmle_MSM$bootstrap.res$beta.A1A2)

    ltmle_MSM$bootstrap.res$RD.A1.A2_0 <- ltmle_MSM$bootstrap.res$p.A1_1.A2_0 - ltmle_MSM$bootstrap.res$p.A1_0.A2_0
    ltmle_MSM$bootstrap.res$RD.A1.A2_1 <- ltmle_MSM$bootstrap.res$p.A1_1.A2_1 - ltmle_MSM$bootstrap.res$p.A1_0.A2_1
    ltmle_MSM$bootstrap.res$RD.A2.A1_0 <- ltmle_MSM$bootstrap.res$p.A1_0.A2_1 - ltmle_MSM$bootstrap.res$p.A1_0.A2_0
    ltmle_MSM$bootstrap.res$RD.A2.A1_1 <- ltmle_MSM$bootstrap.res$p.A1_1.A2_1 - ltmle_MSM$bootstrap.res$p.A1_1.A2_0

    ltmle_MSM$bootstrap.res$lnRR.A1.A2_0 <- log(ltmle_MSM$bootstrap.res$p.A1_1.A2_0) - log(ltmle_MSM$bootstrap.res$p.A1_0.A2_0)
    ltmle_MSM$bootstrap.res$lnRR.A1.A2_1 <- log(ltmle_MSM$bootstrap.res$p.A1_1.A2_1) - log(ltmle_MSM$bootstrap.res$p.A1_0.A2_1)
    ltmle_MSM$bootstrap.res$lnRR.A2.A1_0 <- log(ltmle_MSM$bootstrap.res$p.A1_0.A2_1) - log(ltmle_MSM$bootstrap.res$p.A1_0.A2_0)
    ltmle_MSM$bootstrap.res$lnRR.A2.A1_1 <- log(ltmle_MSM$bootstrap.res$p.A1_1.A2_1) - log(ltmle_MSM$bootstrap.res$p.A1_1.A2_0)

    ltmle_MSM$bootstrap.res$a.INT <- ltmle_MSM$bootstrap.res$p.A1_1.A2_1 -
      ltmle_MSM$bootstrap.res$p.A1_1.A2_0 -
      ltmle_MSM$bootstrap.res$p.A1_0.A2_1 +
      ltmle_MSM$bootstrap.res$p.A1_0.A2_0

    ltmle_MSM$bootstrap.res$RERI <- (ltmle_MSM$bootstrap.res$p.A1_1.A2_1 -
                                       ltmle_MSM$bootstrap.res$p.A1_1.A2_0 -
                                       ltmle_MSM$bootstrap.res$p.A1_0.A2_1 +
                                       ltmle_MSM$bootstrap.res$p.A1_0.A2_0) / ltmle_MSM$bootstrap.res$p.A1_0.A2_0

    ltmle_MSM$bootstrap.res$ln.m.INT <- log(ltmle_MSM$bootstrap.res$p.A1_1.A2_1) + log(ltmle_MSM$bootstrap.res$p.A1_0.A2_0) -
                                              log(ltmle_MSM$bootstrap.res$p.A1_1.A2_0) - log(ltmle_MSM$bootstrap.res$p.A1_0.A2_1)

    # A1 = 0 et A2 = 0
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0] <- sd(ltmle_MSM$bootstrap.res$p.A1_0.A2_0)
    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0]
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0]

    # A1 = 1 et A2 = 0
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0] <- sd(ltmle_MSM$bootstrap.res$p.A1_1.A2_0)
    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0]
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0]

    # A1 = 0 et A2 = 1
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$p.A1_0.A2_1)
    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1]
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 0 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1]

    # A1 = 1 et A2 = 1
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$p.A1_1.A2_1)
    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$p[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1]

    # RD.A1.A2is0
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- sd(ltmle_MSM$bootstrap.res$RD.A1.A2_0)
    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] -
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0]
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] +
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0]

    # RD.A1.A2is1
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$RD.A1.A2_1)
    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1]

    # RD.A2.A1is0
    int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$RD.A2.A1_0)
    int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1]
    int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1]

    # RD.A2.A1is1
    int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$RD.A2.A1_1)
    int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1]

    # RR.A1.A2is0
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- sd(ltmle_MSM$bootstrap.res$lnRR.A1.A2_0)
    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0])
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0])

    # RR.A1.A2is1
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$lnRR.A1.A2_1)
    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                               qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                               qnorm(0.975) * int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1])

    # RR.A2.A1is0
    int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$lnRR.A2.A1_0)
    int.r$RR.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1])
    int.r$RR.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1])

    # RR.A2.A1is1
    int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$lnRR.A2.A1_1)
    int.r$RR.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$RR.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1])

    # additive interaction
    int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$a.INT)
    int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] -
      qnorm(0.975) * int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1]
    int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] +
      qnorm(0.975) * int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1]

    # RERI
    int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$RERI)
    int.r$RERI.lo[int.r$A1 == 1 & int.r$A2 == 1] <- (int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] -
                                                       qnorm(0.975) * int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$RERI.up[int.r$A1 == 1 & int.r$A2 == 1] <- (int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] +
                                                       qnorm(0.975) * int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1])

    # multiplicative interaction
    int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- sd(ltmle_MSM$bootstrap.res$ln.m.INT)
    int.r$m.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1]) -
                                                           qnorm(0.975) * int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1])
    int.r$m.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- exp(log(int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1]) +
                                                           qnorm(0.975) * int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1])

    bootstrap.res <- ltmle_MSM$bootstrap.res
  }


  # back transformation for continous outcomes
  if (ltmle_MSM$ltmle_MSM$transformOutcome == TRUE) {
    back.trfs <- function(y, range) {
      return((y * (range[2] - range[1])) + range[1])
    }

    range.Y <- attr(ltmle_MSM$ltmle_MSM$transformOutcome, "Yrange")

    # A1 = 0 et A2 = 0
    int.r$p[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # A1 = 1 et A2 = 0
    int.r$p[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # A1 = 0 et A2 = 1
    int.r$p[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # A1 = 1 et A2 = 1
    int.r$p[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p[int.r$A1 == 1 & int.r$A2 == 1], range.Y)

    # RD.A1.A2is0
    int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # RD.A1.A2is1
    int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RD.A2.A1is0
    int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RD.A2.A1is1
    int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])

    # RR.A1.A2is0
    int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # RR.A1.A2is1
    int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # RR.A2.A1is0
    int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # RR.A2.A1is1
    int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- NA


    # additive interaction
    int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RERI
    int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # multiplicative interaction
    int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- NA

    # A1 = 0 et A2 = 0
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # A1 = 1 et A2 = 0
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # A1 = 0 et A2 = 1
    int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # A1 = 1 et A2 = 1
    int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1], range.Y)
    int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1], range.Y)
    int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1], range.Y)

    # RD.A1.A2is0
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # RD.A1.A2is1
    int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RD.A2.A1is0
    int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RD.A2.A1is1
    int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])

    # RR.A1.A2is0
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # RR.A1.A2is1
    int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # RR.A2.A1is0
    int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    int.r$RR.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    int.r$RR.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # RR.A2.A1is1
    int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RR.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RR.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA

    # additive interaction
    int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # RERI
    int.r$sd.RERI[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RERI.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$RERI.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # multiplicative interaction
    int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$m.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    int.r$m.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
  }

  return(list(int.r = int.r,
              Anodes = ltmle_MSM$Anodes,
              Ynodes = ltmle_MSM$Ynodes,
              transformOutcome = ltmle_MSM$ltmle_MSM$transformOutcome,
              bootstrap.res = bootstrap.res))
}
