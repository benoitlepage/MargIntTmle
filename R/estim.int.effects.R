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
#' df <- generate.data.multcat(N = 1000, b = param.causal.model.multcat())
#' head(df)
#' df <- data.frame(df[,c("conf1","conf2","conf3")],
#'                  behav.2 = ifelse(df$behav == 2, 1, 0),
#'                  behav.3 = ifelse(df$behav == 3, 1, 0),
#'                  env.2 = ifelse(df$env == 2, 1, 0),
#'                  env.3 = ifelse(df$env == 3, 1, 0),
#'                  hlth.outcome = df$hlth.outcome)
#' head(df)
#'
#' # Define Q and g formulas
#' # an A1 * A2 interaction term is recommended in the Q formula for the estimation
#' # of interaction effects
#' Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
#' g_formulas = c("behav.2 ~ conf1 + conf2",
#'                "behav.3 ~ conf1 + conf2 + behav.2",
#'                "env.2 ~ conf1 + conf3",
#'                "env.3 ~ conf1 + conf3 + env.2")
#'
#' # Define SuperLearner libraries
#' SL.library = list(Q=list("SL.glm", c("SL.glm", "screen.corP"),"SL.glmnet", "SL.mean"),
#'                   g=list("SL.glm", c("SL.glm", "screen.corP"),"SL.glmnet", "SL.mean"))
#'
#' # Estimate MSM parameters by IPTW and TMLE
#' interaction.ltmle <- int.ltmleMSM(data = df,
#'                                   Qform = Q_formulas,
#'                                   gform = g_formulas,
#'                                   A1nodes = c("behav.2","behav.3"),
#'                                   A2nodes = c("env.2","env.3"),
#'                                   Vnodes = NULL,
#'                                   Lnodes = NULL,
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
                              estimator = c("tmle", "iptw", "gcomp")) {

  data <- ltmle_MSM$data

  if(estimator == "gcomp") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == FALSE) stop("The ltmle function did not use the gcomp estimator, but the iptw +/- tmle estimator"))

    beta <- ltmle_MSM$ltmle_MSM$beta
  }

  if(estimator == "iptw") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == TRUE) stop("The ltmle function used the gcomp estimator"))

    beta <- ltmle_MSM$ltmle_MSM$beta.iptw
    IC <- ltmle_MSM$ltmle_MSM$IC.iptw
    var_IC <- var(IC) / nrow(data)
  }

  if(estimator == "tmle") {
    try(if(ltmle_MSM$ltmle_MSM$gcomp == TRUE) stop("The ltmle function used the gcomp estimator, tmle is not available"))

    beta <- ltmle_MSM$ltmle_MSM$beta
    IC <- ltmle_MSM$ltmle_MSM$IC
    var_IC <- var(IC) / nrow(data)
  }

  ## save results in 5 tables: probabilities, risk differences, relative risks, odds ratios and interactions
  # interactions A1 * A2
  if(is.null(ltmle_MSM$Vnodes)) {
    probs <- matrix(NA,
                    nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$A2nodes) + 1),
                    ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$A2nodes) +
                              4)) # p, sd.p, p.lo, p.up
    colnames(probs) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes, "p", "sd.p", "p.lo", "p.up")

    RD <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$A2nodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$A2nodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # RD.A1
                           length(ltmle_MSM$A2nodes) * 4 )) # RD.A2
    RD.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      RD.A1.names <- c(RD.A1.names,
                       paste0("RD.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.RD.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.RD.",ltmle_MSM$A1nodes)[i],
                       paste0("up.RD.",ltmle_MSM$A1nodes)[i])
    }
    RD.A2.names <- NULL
    for(i in 1:length(ltmle_MSM$A2nodes)) {
      RD.A2.names <- c(RD.A2.names,
                       paste0("RD.",ltmle_MSM$A2nodes)[i],
                       paste0("sd.RD.",ltmle_MSM$A2nodes)[i],
                       paste0("lo.RD.",ltmle_MSM$A2nodes)[i],
                       paste0("up.RD.",ltmle_MSM$A2nodes)[i])
    }
    colnames(RD) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                      RD.A1.names,RD.A2.names)
    rm(RD.A1.names,RD.A2.names)


    RR <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$A2nodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$A2nodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # RR.A1
                           length(ltmle_MSM$A2nodes) * 4 )) # RR.A2
    RR.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      RR.A1.names <- c(RR.A1.names,
                       paste0("RR.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.lnRR.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.RR.",ltmle_MSM$A1nodes)[i],
                       paste0("up.RR.",ltmle_MSM$A1nodes)[i])
    }
    RR.A2.names <- NULL
    for(i in 1:length(ltmle_MSM$A2nodes)) {
      RR.A2.names <- c(RR.A2.names,
                       paste0("RR.",ltmle_MSM$A2nodes)[i],
                       paste0("sd.lnRR.",ltmle_MSM$A2nodes)[i],
                       paste0("lo.RR.",ltmle_MSM$A2nodes)[i],
                       paste0("up.RR.",ltmle_MSM$A2nodes)[i])
    }
    colnames(RR) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                      RR.A1.names,RR.A2.names)
    rm(RR.A1.names,RR.A2.names)

    OR <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$A2nodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$A2nodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # OR.A1
                           length(ltmle_MSM$A2nodes) * 4 )) # OR.A2
    OR.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      OR.A1.names <- c(OR.A1.names,
                       paste0("OR.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.lnOR.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.OR.",ltmle_MSM$A1nodes)[i],
                       paste0("up.OR.",ltmle_MSM$A1nodes)[i])
    }
    OR.A2.names <- NULL
    for(i in 1:length(ltmle_MSM$A2nodes)) {
      OR.A2.names <- c(OR.A2.names,
                       paste0("OR.",ltmle_MSM$A2nodes)[i],
                       paste0("sd.lnOR.",ltmle_MSM$A2nodes)[i],
                       paste0("lo.OR.",ltmle_MSM$A2nodes)[i],
                       paste0("up.OR.",ltmle_MSM$A2nodes)[i])
    }
    colnames(OR) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                      OR.A1.names,OR.A2.names)
    rm(OR.A1.names,OR.A2.names)

    int <- matrix(NA,
                  nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$A2nodes) + 1),
                  ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$A2nodes) +
                            3 * 4) ) # a.INT, RERI, m.INT
    colnames(int) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                       "a.INT", "sd.a.INT", "lo.a.INT", "up.a.INT",
                       "RERI","sd.RERI","lo.RERI","up.RERI",
                       "m.INT", "sd.ln.m.INT", "lo.m.INT", "up.m.INT")
  }

  # effect of A1, modified by V
  if(is.null(ltmle_MSM$A2nodes)) {
    probs <- matrix(NA,
                    nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$Vnodes) + 1),
                    ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$Vnodes) +
                              4)) # p, sd.p, p.lo, p.up
    colnames(probs) <- c(ltmle_MSM$A1nodes, ltmle_MSM$Vnodes, "p", "sd.p", "p.lo", "p.up")

    RD <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$Vnodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$Vnodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # RD.A1
                           length(ltmle_MSM$Vnodes) * 4 )) # RD.V
    RD.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      RD.A1.names <- c(RD.A1.names,
                       paste0("RD.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.RD.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.RD.",ltmle_MSM$A1nodes)[i],
                       paste0("up.RD.",ltmle_MSM$A1nodes)[i])
    }
    RD.V.names <- NULL
    for(i in 1:length(ltmle_MSM$Vnodes)) {
      RD.V.names <- c(RD.V.names,
                       paste0("RD.",ltmle_MSM$Vnodes)[i],
                       paste0("sd.RD.",ltmle_MSM$Vnodes)[i],
                       paste0("lo.RD.",ltmle_MSM$Vnodes)[i],
                       paste0("up.RD.",ltmle_MSM$Vnodes)[i])
    }
    colnames(RD) <- c(ltmle_MSM$A1nodes, ltmle_MSM$Vnodes,
                      RD.A1.names,RD.V.names)
    rm(RD.A1.names,RD.V.names)


    RR <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$Vnodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$Vnodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # RR.A1
                           length(ltmle_MSM$Vnodes) * 4 )) # RR.A2
    RR.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      RR.A1.names <- c(RR.A1.names,
                       paste0("RR.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.lnRR.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.RR.",ltmle_MSM$A1nodes)[i],
                       paste0("up.RR.",ltmle_MSM$A1nodes)[i])
    }
    RR.V.names <- NULL
    for(i in 1:length(ltmle_MSM$Vnodes)) {
      RR.V.names <- c(RR.V.names,
                       paste0("RR.",ltmle_MSM$Vnodes)[i],
                       paste0("sd.lnRR.",ltmle_MSM$Vnodes)[i],
                       paste0("lo.RR.",ltmle_MSM$Vnodes)[i],
                       paste0("up.RR.",ltmle_MSM$Vnodes)[i])
    }
    colnames(RR) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                      RR.A1.names,RR.V.names)
    rm(RR.A1.names,RR.V.names)

    OR <- matrix(NA,
                 nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$Vnodes) + 1),
                 ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$Vnodes) +
                           length(ltmle_MSM$A1nodes) * 4 + # OR.A1
                           length(ltmle_MSM$Vnodes) * 4 )) # OR.A2
    OR.A1.names <- NULL
    for(i in 1:length(ltmle_MSM$A1nodes)) {
      OR.A1.names <- c(OR.A1.names,
                       paste0("OR.",ltmle_MSM$A1nodes)[i],
                       paste0("sd.lnOR.",ltmle_MSM$A1nodes)[i],
                       paste0("lo.OR.",ltmle_MSM$A1nodes)[i],
                       paste0("up.OR.",ltmle_MSM$A1nodes)[i])
    }
    OR.V.names <- NULL
    for(i in 1:length(ltmle_MSM$Vnodes)) {
      OR.V.names <- c(OR.V.names,
                      paste0("OR.",ltmle_MSM$Vnodes)[i],
                      paste0("sd.lnOR.",ltmle_MSM$Vnodes)[i],
                      paste0("lo.OR.",ltmle_MSM$Vnodes)[i],
                      paste0("up.OR.",ltmle_MSM$Vnodes)[i])
    }
    colnames(OR) <- c(ltmle_MSM$A1nodes, ltmle_MSM$A2nodes,
                      OR.A1.names,OR.V.names)
    rm(OR.A1.names,OR.V.names)

    int <- matrix(NA,
                  nrow = (length(ltmle_MSM$A1nodes) + 1) * (length(ltmle_MSM$Vnodes) + 1),
                  ncol = (length(ltmle_MSM$A1nodes) + length(ltmle_MSM$Vnodes) +
                            3 * 4) ) # a.INT, RERI, m.INT
    colnames(int) <- c(ltmle_MSM$A1nodes, ltmle_MSM$Vnodes,
                       "a.INT", "sd.a.INT", "lo.a.INT", "up.a.INT",
                       "RERI","sd.RERI","lo.RERI","up.RERI",
                       "m.INT", "sd.ln.m.INT", "lo.m.INT", "up.m.INT")
  }

  values.A1 <- matrix(0, nrow = length(ltmle_MSM$A1nodes) + 1, ncol = length(ltmle_MSM$A1nodes))
  for(i in 1:length(ltmle_MSM$A1nodes)) {
    values.A1[i + 1, i] <- 1
    colnames(values.A1) <- ltmle_MSM$A1nodes
  }

  if(is.null(ltmle_MSM$Vnodes)) {
    values.A2 <- matrix(0, nrow = length(ltmle_MSM$A2nodes) + 1, ncol = length(ltmle_MSM$A2nodes))
    for(i in 1:length(ltmle_MSM$A2nodes)) {
      values.A2[i + 1, i] <- 1
      colnames(values.A2) <- ltmle_MSM$A2nodes
    }
  }
  if(is.null(ltmle_MSM$A2nodes)) {
    values.A2 <- matrix(0, nrow = length(ltmle_MSM$Vnodes) + 1, ncol = length(ltmle_MSM$Vnodes))
    for(i in 1:length(ltmle_MSM$Vnodes)) {
      values.A2[i + 1, i] <- 1
      colnames(values.A2) <- ltmle_MSM$Vnodes
    }
  }

  row <- 1
  for(i in 1:nrow(values.A1)) {
    for(j in 1:nrow(values.A2)) {
      probs[row,1:(ncol(values.A1) + ncol(values.A2))] <- c(values.A1[i,],values.A2[j,])
      RD[row,1:(ncol(values.A1) + ncol(values.A2))] <- c(values.A1[i,],values.A2[j,])
      RR[row,1:(ncol(values.A1) + ncol(values.A2))] <- c(values.A1[i,],values.A2[j,])
      int[row,1:(ncol(values.A1) + ncol(values.A2))] <- c(values.A1[i,],values.A2[j,])

      row <- row + 1
    }
  }
  regimes <- probs[,1:(ncol(values.A1) + ncol(values.A2))]
  Xtemp <- model.matrix(as.formula(ltmle_MSM$working.msm), data = data.frame(regimes, Y = 1))

  if(estimator == "tmle" | estimator = "iptw") {
    ## probabilities
    probs[,"p"] <- plogis(Xtemp %*% beta)

    # grad.p <- list()
    for(i in 1:nrow(probs)) {
      grad.p <- Xtemp[i,] * probs[i,"p"] * (1 - probs[i,"p"])
      probs[i,"sd.p"] <- sqrt(t(grad.p) %*% var_IC %*% grad.p)
      probs[i,"p.lo"] <- probs[i,"p"] - qnorm(0.975) * probs[i,"sd.p"]
      probs[i,"p.up"] <- probs[i,"p"] + qnorm(0.975) * probs[i,"sd.p"]
    }


    ## risk differences RD
    # RD.A1|A2
    for(i in 1:ncol(values.A1)) {
      # column index
      c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A2)) {
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
        RD[r.ind, c.RD] <- probs[r.ind,"p"] - probs[r.ref,"p"]
        grad.RD <- (Xtemp[r.ind,] * probs[r.ind,"p"] * (1 - probs[r.ind,"p"]) -
                      Xtemp[r.ref,] * probs[r.ref,"p"] * (1 - probs[r.ref,"p"]))
        RD[r.ind, c.RD + 1] <- sqrt(t(grad.RD) %*% var_IC %*% grad.RD) # sd
        RD[r.ind, c.RD + 2] <- RD[r.ind, c.RD] - qnorm(0.975) * RD[r.ind, c.RD + 1]
        RD[r.ind, c.RD + 3] <- RD[r.ind, c.RD] + qnorm(0.975) * RD[r.ind, c.RD + 1]
      }
    }
    # RD.A2|A1
    for(i in 1:ncol(values.A2)) {
      # column index
      c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A1)) {
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
        RD[r.ind, c.RD] <- probs[r.ind,"p"] - probs[r.ref,"p"]
        grad.RD <- (Xtemp[r.ind,] * probs[r.ind,"p"] * (1 - probs[r.ind,"p"]) -
                      Xtemp[r.ref,] * probs[r.ref,"p"] * (1 - probs[r.ref,"p"]))
        RD[r.ind, c.RD + 1] <- sqrt(t(grad.RD) %*% var_IC %*% grad.RD) # sd
        RD[r.ind, c.RD + 2] <- RD[r.ind, c.RD] - qnorm(0.975) * RD[r.ind, c.RD + 1]
        RD[r.ind, c.RD + 3] <- RD[r.ind, c.RD] + qnorm(0.975) * RD[r.ind, c.RD + 1]
      }
    }

    ## relative risk RR
    # RR.A1|A2
    for(i in 1:ncol(values.A1)) {
      # column index
      c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A2)) {
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
        RR[r.ind, c.RR] <- probs[r.ind,"p"] / probs[r.ref,"p"]
        lnRR <- log(probs[r.ind,"p"]) - log(probs[r.ref,"p"])
        grad.lnRR <- (Xtemp[r.ind,] * (1 - probs[r.ind,"p"]) -
                        Xtemp[r.ref,] * (1 - probs[r.ref,"p"]))
        RR[r.ind, c.RR + 1] <- sqrt(t(grad.lnRR) %*% var_IC %*% grad.lnRR) # sd
        RR[r.ind, c.RR + 2] <- exp(lnRR - qnorm(0.975) * RR[r.ind, c.RR + 1])
        RR[r.ind, c.RR + 3] <- exp(lnRR + qnorm(0.975) * RR[r.ind, c.RR + 1])
      }
    }
    # RR.A2|A1
    for(i in 1:ncol(values.A2)) {
      # column index
      c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A1)) {
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
        RR[r.ind, c.RR] <- probs[r.ind,"p"] / probs[r.ref,"p"]
        lnRR <- log(probs[r.ind,"p"]) - log(probs[r.ref,"p"])
        grad.lnRR <- (Xtemp[r.ind,] * (1 - probs[r.ind,"p"]) -
                        Xtemp[r.ref,] * (1 - probs[r.ref,"p"]))
        RR[r.ind, c.RR + 1] <- sqrt(t(grad.lnRR) %*% var_IC %*% grad.lnRR) # sd
        RR[r.ind, c.RR + 2] <- exp(lnRR - qnorm(0.975) * RR[r.ind, c.RR + 1])
        RR[r.ind, c.RR + 3] <- exp(lnRR + qnorm(0.975) * RR[r.ind, c.RR + 1])
      }
    }

    ## odds ratio
    # OR.A1|A2
    for(i in 1:ncol(values.A1)) {
      # column index
      c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A2)) {
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
        OR[r.ind, c.OR] <- (probs[r.ind,"p"] / (1 - probs[r.ind,"p"])) / (probs[r.ref,"p"] / (1 - probs[r.ref,"p"]))
        lnOR <- log(OR[r.ind, c.OR])
        grad.lnOR <- Xtemp[r.ind,] - Xtemp[r.ref,]
        OR[r.ind, c.OR + 1] <- sqrt(t(grad.lnOR) %*% var_IC %*% grad.lnOR) # sd
        OR[r.ind, c.OR + 2] <- exp(lnOR - qnorm(0.975) * OR[r.ind, c.OR + 1])
        OR[r.ind, c.OR + 3] <- exp(lnOR + qnorm(0.975) * OR[r.ind, c.OR + 1])
      }
    }
    # OR.A2|A1
    for(i in 1:ncol(values.A2)) {
      # column index
      c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
      for(j in 1:nrow(values.A1)) {
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
          cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
        }
        if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
          cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
        }
        # row of the index probability
        r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
        # row of the reference probability
        r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
        OR[r.ind, c.OR] <- (probs[r.ind,"p"] / (1 - probs[r.ind,"p"])) / (probs[r.ref,"p"] / (1 - probs[r.ref,"p"]))
        lnOR <- log(OR[r.ind, c.OR])
        grad.lnOR <- Xtemp[r.ind,] - Xtemp[r.ref,]
        OR[r.ind, c.OR + 1] <- sqrt(t(grad.lnOR) %*% var_IC %*% grad.lnOR) # sd
        OR[r.ind, c.OR + 2] <- exp(lnOR - qnorm(0.975) * OR[r.ind, c.OR + 1])
        OR[r.ind, c.OR + 3] <- exp(lnOR + qnorm(0.975) * OR[r.ind, c.OR + 1])
      }
    }

    ## interactions
    for(i in 1:ncol(values.A1)) {
      for(j in 1:ncol(values.A2)) {
        # define relevant row indices in probs table
        r.00 <- which(rowSums(regimes[,c(colnames(values.A1),colnames(values.A2))]) == 0)
        r.10 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                        (rowSums(as.matrix(regimes[,c(colnames(values.A2))])) == 0))
        r.01 <- which((rowSums(as.matrix(regimes[,c(colnames(values.A1))])) == 0) &
                        (regimes[,c(colnames(values.A2)[j])] == 1))
        r.11 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                        (regimes[,c(colnames(values.A2)[j])] == 1))

        int[r.11,"a.INT"] <- probs[r.11,"p"] - probs[r.10,"p"] - probs[r.01,"p"] + probs[r.00,"p"]
        grad.a.INT <- (Xtemp[r.11,] * probs[r.11,"p"] * (1 - probs[r.11,"p"]) -
                         Xtemp[r.10,] * probs[r.10,"p"] * (1 - probs[r.10,"p"]) -
                         Xtemp[r.01,] * probs[r.01,"p"] * (1 - probs[r.01,"p"]) +
                         Xtemp[r.00,] * probs[r.00,"p"] * (1 - probs[r.00,"p"]))
        int[r.11, "sd.a.INT"] <- sqrt(t(grad.a.INT) %*% var_IC %*% grad.a.INT) # sd
        int[r.11, "lo.a.INT"] <- int[r.11,"a.INT"] - qnorm(0.975) * int[r.11, "sd.a.INT"]
        int[r.11, "up.a.INT"] <- int[r.11,"a.INT"] + qnorm(0.975) * int[r.11, "sd.a.INT"]

        int[r.11,"RERI"] <- (probs[r.11,"p"] - probs[r.10,"p"] - probs[r.01,"p"] + probs[r.00,"p"]) / probs[r.00,"p"]
        grad.RERI <- (grad.a.INT / probs[r.00,"p"]) -
          (probs[r.11,"p"] - probs[r.10,"p"] - probs[r.01,"p"] + probs[r.00,"p"]) * (1 - probs[r.00,"p"]) / probs[r.00,"p"]
        int[r.11, "sd.RERI"] <- sqrt(t(grad.RERI) %*% var_IC %*% grad.RERI) # sd
        int[r.11, "lo.RERI"] <- int[r.11,"RERI"] - qnorm(0.975) * int[r.11, "sd.RERI"]
        int[r.11, "up.RERI"] <- int[r.11,"RERI"] + qnorm(0.975) * int[r.11, "sd.RERI"]

        int[r.11,"m.INT"] <- (probs[r.11,"p"] * probs[r.00,"p"]) / (probs[r.10,"p"] * probs[r.01,"p"])
        ln.m.INT <- log(probs[r.11,"p"]) - log(probs[r.10,"p"]) - log(probs[r.01,"p"]) + log(probs[r.00,"p"])
        grad.ln.m.INT <- (Xtemp[r.11,] * (1 - probs[r.11,"p"]) -
                            Xtemp[r.10,] * (1 - probs[r.10,"p"]) -
                            Xtemp[r.01,] * (1 - probs[r.01,"p"]) +
                            Xtemp[r.00,] * (1 - probs[r.00,"p"]))
        int[r.11, "sd.ln.m.INT"] <- sqrt(t(grad.ln.m.INT) %*% var_IC %*% grad.ln.m.INT) # sd
        int[r.11, "lo.m.INT"] <- exp(ln.m.INT - qnorm(0.975) * int[r.11, "sd.ln.m.INT"])
        int[r.11, "up.m.INT"] <- exp(ln.m.INT + qnorm(0.975) * int[r.11, "sd.ln.m.INT"])
      }
    }

    ## bootstrap
    bootstrap.res <- ltmle_MSM$bootstrap.res
  }


  if(estimator == "gcomp") {
    B <- nrow(ltmle_MSM$bootstrap.res)

    ## probabilities
    probs[,"p"] <- plogis(Xtemp %*% beta)
    probs.boot <- array(NA, dim = c(dim(probs[,c(colnames(values.A1),colnames(values.A2),"p")]),B))
    colnames(probs.boot) <- c(colnames(values.A1),colnames(values.A2),"p")
    probs.boot[,c(colnames(values.A1),colnames(values.A2)),] <- regimes
    for(b in 1:B) {
      beta.boot <- ltmle_MSM$bootstrap.res[b,]
      probs.boot[,"p",b] <- plogis(Xtemp %*% as.numeric(ltmle_MSM$bootstrap.res[b,]))
    }
    for(r in 1:nrow(probs)) {
      probs[r,"sd.p"] <- sd(probs.boot[r,"p",])
      probs[r,"p.lo"] <- probs[r,"p"] - qnorm(0.975) * sd(probs.boot[r,"p",])
      probs[r,"p.up"] <- probs[r,"p"] + qnorm(0.975) * sd(probs.boot[r,"p",])
    }


    ## Risk difference
    RD.boot <- array(NA, dim = c(dim(RD),B))
    colnames(RD.boot) <- colnames(RD)
    RD.boot[,c(colnames(values.A1),colnames(values.A2)),] <- regimes
    for(b in 1:B) {
      # RD.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          RD.boot[r.ind, c.RD, b] <- probs.boot[r.ind,"p",b] - probs.boot[r.ref,"p",b]
        }
      }
      # RD.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          RD.boot[r.ind, c.RD,b] <- probs.boot[r.ind,"p",b] - probs.boot[r.ref,"p",b]
        }
      }
    }

    for(r in 1:nrow(RD)) {
      # RD.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          RD[r.ind, c.RD] <- probs[r.ind,"p"] - probs[r.ref,"p"]
          RD[r.ind, c.RD + 1] <- sd(RD.boot[r.ind, c.RD,])
          RD[r.ind, c.RD + 2] <- RD[r.ind, c.RD] - qnorm(0.975) * RD[r.ind, c.RD + 1]
          RD[r.ind, c.RD + 3] <- RD[r.ind, c.RD] + qnorm(0.975) * RD[r.ind, c.RD + 1]
        }
      }
      # RD.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.RD <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          RD[r.ind, c.RD] <- probs[r.ind,"p"] - probs[r.ref,"p"]
          RD[r.ind, c.RD + 1] <- sd(RD.boot[r.ind, c.RD,])
          RD[r.ind, c.RD + 2] <- RD[r.ind, c.RD] - qnorm(0.975) * RD[r.ind, c.RD + 1]
          RD[r.ind, c.RD + 3] <- RD[r.ind, c.RD] + qnorm(0.975) * RD[r.ind, c.RD + 1]
        }
      }
    }

    ## Relative risk
    RR.boot <- array(NA, dim = c(dim(RR),B))
    colnames(RR.boot) <- colnames(RR)
    RR.boot[,c(colnames(values.A1),colnames(values.A2)),] <- regimes
    for(b in 1:B) {
      # RR.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          RR.boot[r.ind, c.RR, b] <- probs.boot[r.ind,"p",b] / probs.boot[r.ref,"p",b]
        }
      }
      # RR.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          RR.boot[r.ind, c.RR,b] <- probs.boot[r.ind,"p",b] / probs.boot[r.ref,"p",b]
        }
      }
    }

    for(r in 1:nrow(RR)) {
      # RR.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          RR[r.ind, c.RR] <- probs[r.ind,"p"] / probs[r.ref,"p"]
          RR[r.ind, c.RR + 1] <- sd(log(RR.boot[r.ind, c.RR,]))
          RR[r.ind, c.RR + 2] <- exp(log(RR[r.ind, c.RR]) - qnorm(0.975) * RR[r.ind, c.RR + 1])
          RR[r.ind, c.RR + 3] <- exp(log(RR[r.ind, c.RR]) + qnorm(0.975) * RR[r.ind, c.RR + 1])
        }
      }
      # RR.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.RR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          RR[r.ind, c.RR] <- probs[r.ind,"p"] / probs[r.ref,"p"]
          RR[r.ind, c.RR + 1] <- sd(log(RR.boot[r.ind, c.RR,]))
          RR[r.ind, c.RR + 2] <- exp(log(RR[r.ind, c.RR]) - qnorm(0.975) * RR[r.ind, c.RR + 1])
          RR[r.ind, c.RR + 3] <- exp(log(RR[r.ind, c.RR]) + qnorm(0.975) * RR[r.ind, c.RR + 1])
        }
      }
    }

    ## Odds ratios
    OR.boot <- array(NA, dim = c(dim(OR),B))
    colnames(OR.boot) <- colnames(OR)
    OR.boot[,c(colnames(values.A1),colnames(values.A2)),] <- regimes
    for(b in 1:B) {
      # OR.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          OR.boot[r.ind, c.OR, b] <- (probs.boot[r.ind,"p",b] / (1 - probs.boot[r.ind,"p",b])) / (probs.boot[r.ref,"p",b] / (1 - probs.boot[r.ref,"p",b]))
        }
      }
      # OR.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          OR.boot[r.ind, c.OR,b] <- (probs.boot[r.ind,"p",b] / (1 - probs.boot[r.ind,"p",b])) / (probs.boot[r.ref,"p",b] / (1 - probs.boot[r.ref,"p",b]))
        }
      }
    }

    for(r in 1:nrow(OR)) {
      # OR.A1|A2
      for(i in 1:ncol(values.A1)) {
        # column index
        c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A2)) {
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A2)]) == values.A2[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A2)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A2)]), 1, function(r) r == values.A2[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A1)])) == 0) & cond2)
          OR[r.ind, c.OR] <- (probs[r.ind,"p"] / (1 - probs[r.ind,"p"])) / (probs[r.ref,"p"] / (1 - probs[r.ref,"p"]))
          OR[r.ind, c.OR + 1] <- sd(log(OR.boot[r.ind, c.OR,]))
          OR[r.ind, c.OR + 2] <- exp(log(OR[r.ind, c.OR]) - qnorm(0.975) * OR[r.ind, c.OR + 1])
          OR[r.ind, c.OR + 3] <- exp(log(OR[r.ind, c.OR]) + qnorm(0.975) * OR[r.ind, c.OR + 1])
        }
      }
      # OR.A2|A1
      for(i in 1:ncol(values.A2)) {
        # column index
        c.OR <- ncol(values.A1) + ncol(values.A2) + 4 * ncol(values.A1) + 4 * (i - 1) + 1
        for(j in 1:nrow(values.A1)) {
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] == 1) {
            cond2 = as.matrix(regimes[,colnames(values.A1)]) == values.A1[j,]
          }
          if(dim(as.matrix(regimes[,colnames(values.A1)]))[2] > 1) {
            cond2 = as.matrix(apply(apply(as.matrix(regimes[,colnames(values.A1)]), 1, function(r) r == values.A1[j,]),2,prod) == 1)
          }
          # row of the index probability
          r.ind <- which(as.matrix(regimes[,ncol(values.A2) + i] == 1) & cond2)
          # row of the reference probability
          r.ref <- which(as.matrix(rowSums(as.matrix(regimes[,colnames(values.A2)])) == 0) & cond2)
          OR[r.ind, c.OR] <- (probs[r.ind,"p"] / (1 - probs[r.ind,"p"])) / (probs[r.ref,"p"] / (1 - probs[r.ref,"p"]))
          OR[r.ind, c.OR + 1] <- sd(log(OR.boot[r.ind, c.OR,]))
          OR[r.ind, c.OR + 2] <- exp(log(OR[r.ind, c.OR]) - qnorm(0.975) * OR[r.ind, c.OR + 1])
          OR[r.ind, c.OR + 3] <- exp(log(OR[r.ind, c.OR]) + qnorm(0.975) * OR[r.ind, c.OR + 1])
        }
      }
    }


    ## interactions
    int.boot <- array(NA, dim = c(dim(int),B))
    colnames(int.boot) <- colnames(int)
    int.boot[,c(colnames(values.A1),colnames(values.A2)),] <- regimes
    for(b in 1:B) {
      for(i in 1:ncol(values.A1)) {
        for(j in 1:ncol(values.A2)) {
          # define relevant row indices in probs table
          r.00 <- which(rowSums(regimes[,c(colnames(values.A1),colnames(values.A2))]) == 0)
          r.10 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                          (rowSums(as.matrix(regimes[,c(colnames(values.A2))])) == 0))
          r.01 <- which((rowSums(as.matrix(regimes[,c(colnames(values.A1))])) == 0) &
                          (regimes[,c(colnames(values.A2)[j])] == 1))
          r.11 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                          (regimes[,c(colnames(values.A2)[j])] == 1))

          int.boot[r.11,"a.INT",b] <- probs.boot[r.11,"p",b] - probs.boot[r.10,"p",b] - probs.boot[r.01,"p",b] + probs.boot[r.00,"p",b]

          int.boot[r.11,"RERI",b] <- (probs.boot[r.11,"p",b] - probs.boot[r.10,"p",b] - probs.boot[r.01,"p",b] + probs.boot[r.00,"p",b]) / probs.boot[r.00,"p",b]

          int.boot[r.11,"m.INT",b] <- (probs.boot[r.11,"p",b] * probs.boot[r.00,"p",b]) / (probs.boot[r.10,"p",b] * probs.boot[r.01,"p",b])
        }
      }
    }

    for(i in 1:ncol(values.A1)) {
      for(j in 1:ncol(values.A2)) {
        # define relevant row indices in probs table
        r.00 <- which(rowSums(regimes[,c(colnames(values.A1),colnames(values.A2))]) == 0)
        r.10 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                        (rowSums(as.matrix(regimes[,c(colnames(values.A2))])) == 0))
        r.01 <- which((rowSums(as.matrix(regimes[,c(colnames(values.A1))])) == 0) &
                        (regimes[,c(colnames(values.A2)[j])] == 1))
        r.11 <- which((regimes[,c(colnames(values.A1)[i])] == 1) &
                        (regimes[,c(colnames(values.A2)[j])] == 1))

        int[r.11,"a.INT"] <- probs[r.11,"p"] - probs[r.10,"p"] - probs[r.01,"p"] + probs[r.00,"p"]
        int[r.11, "sd.a.INT"] <- sd(int.boot[r.11,"a.INT",])
        int[r.11, "lo.a.INT"] <- int[r.11,"a.INT"] - qnorm(0.975) * int[r.11, "sd.a.INT"]
        int[r.11, "up.a.INT"] <- int[r.11,"a.INT"] + qnorm(0.975) * int[r.11, "sd.a.INT"]

        int[r.11,"RERI"] <- (probs[r.11,"p"] - probs[r.10,"p"] - probs[r.01,"p"] + probs[r.00,"p"]) / probs[r.00,"p"]
        int[r.11, "sd.RERI"] <- sd(int.boot[r.11,"RERI",])
        int[r.11, "lo.RERI"] <- int[r.11,"RERI"] - qnorm(0.975) * int[r.11, "sd.RERI"]
        int[r.11, "up.RERI"] <- int[r.11,"RERI"] + qnorm(0.975) * int[r.11, "sd.RERI"]

        int[r.11,"m.INT"] <- (probs[r.11,"p"] * probs[r.00,"p"]) / (probs[r.10,"p"] * probs[r.01,"p"])
        ln.m.INT <- log(probs[r.11,"p"]) - log(probs[r.10,"p"]) - log(probs[r.01,"p"]) + log(probs[r.00,"p"])
        int[r.11, "sd.ln.m.INT"] <- sd(log(int.boot[r.11,"m.INT",]))
        int[r.11, "lo.m.INT"] <- exp(ln.m.INT - qnorm(0.975) * int[r.11, "sd.ln.m.INT"])
        int[r.11, "up.m.INT"] <- exp(ln.m.INT + qnorm(0.975) * int[r.11, "sd.ln.m.INT"])
      }
    }

    bootstrap.res <- list(probs.boot = probs.boot,
                          RD.boot = RD.boot[,c(colnames(values.A1),colnames(values.A2),
                                               paste0("RD.",colnames(values.A1)),
                                               paste0("RD.",colnames(values.A2))),],
                          RR.boot = RR.boot[,c(colnames(values.A1),colnames(values.A2),
                                               paste0("RR.",colnames(values.A1)),
                                               paste0("RR.",colnames(values.A2))),],
                          OR.boot = OR.boot[,c(colnames(values.A1),colnames(values.A2),
                                               paste0("OR.",colnames(values.A1)),
                                               paste0("OR.",colnames(values.A2))),],
                          int.boot = int.boot[,c(colnames(values.A1),colnames(values.A2),
                                                 "a.INT","RERI","m.INT"),])
  }


  # back transformation for continous outcomes
  if (ltmle_MSM$ltmle_MSM$transformOutcome == TRUE) {
    back.trfs <- function(y, range) {
      return(ifelse(!is.na(y),
                    (y * (range[2] - range[1])) + range[1],
                    NA))
    }

    range.Y <- attr(ltmle_MSM$ltmle_MSM$transformOutcome, "Yrange")

    for(r in 1:nrow(probs)) {
      for(c in (ncol(values.A1) + ncol(values.A2) + 1):ncol(probs)) {
        probs[r,c] <- back.trfs(probs[r,c], range.Y)
      }
    }

    for(r in 1:nrow(RD)) {
      for(c in (ncol(values.A1) + ncol(values.A2) + 1):ncol(RD)) {
        RD[r,c] <- RD[r,c] * (range.Y[2] - range.Y[1])
      }
    }

    RR[,(ncol(values.A1) + ncol(values.A2) + 1):ncol(RD)] <- NA
    OR[,(ncol(values.A1) + ncol(values.A2) + 1):ncol(RD)] <- NA

    # # A1 = 0 et A2 = 0
    # int.r$p[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # # A1 = 1 et A2 = 0
    # int.r$p[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # # A1 = 0 et A2 = 1
    # int.r$p[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # # A1 = 1 et A2 = 1
    # int.r$p[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p[int.r$A1 == 1 & int.r$A2 == 1], range.Y)

    # # RD.A1.A2is0
    # int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # # RD.A1.A2is1
    # int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RD.A2.A1is0
    # int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RD.A2.A1is1
    # int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])

    # # RR.A1.A2is0
    # int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # # RR.A1.A2is1
    # int.r$RR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # # RR.A2.A1is0
    # int.r$RR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # # RR.A2.A1is1
    # int.r$RR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- NA

    for(r in 1:nrow(probs)) {
      int[r,"a.INT"] <- int[r,"a.INT"] * (range.Y[2] - range.Y[1])
      int[r,"sd.a.INT"] <- int[r,"sd.a.INT"] * (range.Y[2] - range.Y[1])
      int[r,"lo.a.INT"] <- int[r,"lo.a.INT"] * (range.Y[2] - range.Y[1])
      int[r,"up.a.INT"] <- int[r,"up.a.INT"] * (range.Y[2] - range.Y[1])
    }
    int[,c("RERI","sd.RERI","lo.RERI","up.RERI","m.INT","sd.ln.m.INT","lo.m.INT","up.m.INT")] <- NA

    # # additive interaction
    # int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RERI
    # int.r$RERI[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # # multiplicative interaction
    # int.r$m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- NA

    # # A1 = 0 et A2 = 0
    # int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0] <- back.trfs(int.r$p.up[int.r$A1 == 0 & int.r$A2 == 0], range.Y)
    # # A1 = 1 et A2 = 0
    # int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0] <- back.trfs(int.r$p.up[int.r$A1 == 1 & int.r$A2 == 0], range.Y)
    # # A1 = 0 et A2 = 1
    # int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$sd.p[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p.lo[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1] <- back.trfs(int.r$p.up[int.r$A1 == 0 & int.r$A2 == 1], range.Y)
    # # A1 = 1 et A2 = 1
    # int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$sd.p[int.r$A1 == 1 & int.r$A2 == 1], range.Y)
    # int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p.lo[int.r$A1 == 1 & int.r$A2 == 1], range.Y)
    # int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1] <- back.trfs(int.r$p.up[int.r$A1 == 1 & int.r$A2 == 1], range.Y)
    #
    # # RD.A1.A2is0
    # int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 0] * (range.Y[2] - range.Y[1])
    # # RD.A1.A2is1
    # int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.RD.A1[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A1.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RD.A2.A1is0
    # int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$sd.RD.A2[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- int.r$RD.A2.up[int.r$A1 == 0 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RD.A2.A1is1
    # int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.RD.A2[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$RD.A2.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    #
    # # RR.A1.A2is0
    # int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 0] <- NA
    # # RR.A1.A2is1
    # int.r$sd.lnRR.A1[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RR.A1.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RR.A1.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # # RR.A2.A1is0
    # int.r$sd.lnRR.A2[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # int.r$RR.A2.lo[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # int.r$RR.A2.up[int.r$A1 == 0 & int.r$A2 == 1] <- NA
    # # RR.A2.A1is1
    # int.r$sd.lnRR.A2[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RR.A2.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RR.A2.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    #
    # # additive interaction
    # int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$sd.a.INT[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- int.r$a.INT.up[int.r$A1 == 1 & int.r$A2 == 1] * (range.Y[2] - range.Y[1])
    # # RERI
    # int.r$sd.lnRERI[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RERI.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$RERI.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # # multiplicative interaction
    # int.r$sd.ln.m.INT[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$m.INT.lo[int.r$A1 == 1 & int.r$A2 == 1] <- NA
    # int.r$m.INT.up[int.r$A1 == 1 & int.r$A2 == 1] <- NA
  }

  return(list(probs = probs, RD = RD, RR = RR, OR = OR, int = int,
              A1nodes = ltmle_MSM$A1nodes,
              A2nodes = ltmle_MSM$A2nodes,
              Vnodes = ltmle_MSM$Vnodes,
              Ynodes = ltmle_MSM$Ynodes,
              transformOutcome = ltmle_MSM$ltmle_MSM$transformOutcome,
              bootstrap.res = bootstrap.res))
}
