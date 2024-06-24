#' Create table of the marginal interaction effects
#'
#' @param int.r an \code{int.r} data frame from the output obtained by the \code{estim.int.effects} function
#' @param probab.digits integer indicating the number of decimal places to be used for probabilities and risk differences (with \code{round()} function)
#' @param RR.digits integer indicating the number of decimal places to be used for relative risks, RERI and multiplicative interaction effect (with \code{round()} function)
#'
#' @return \code{out.int.table} returns a list of 2 objects:
#'          \itemize{ \item \code{out.table} a data frame with the estimation of various quantities of interest for the interaction effect of A1 * A2 -> Y. The table is presented following the recommendations of Knol and VanderWeele (2012)
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
#' # Estimate quantities of interest for the interaction effect of A1 * A2 -> Y
#' interaction.det <- estim.int.effects(interaction.ltmle, estimator = "tmle")
#'
#' # Show results in a table and give interaction effects
#' table.interactions <- out.int.table(int.r = interaction.det)
#' table.interactions$out.table
#' table.interactions$interaction.effects
out.int.table <- function(int.r = int.r,
                          probab.digits = 3,
                          RR.OR.digits = 2,
                          multipl = c("RR"), # "RR" or "OR"
                          labels.A1 = NULL,
                          labels.A2 = NULL,
                          labels.V = NULL) {
  if(is.null(int.r$A2nodes)) {
    A2nodes <- int.r$Vnodes
    labels.A2 <- labels.V
    values.A2 <- int.r$values.V
  }
  if(is.null(int.r$Vnodes)) {
    A2nodes <- int.r$A2nodes
    values.A2 <- int.r$values.A2
  }

  if(multipl == "RR") {
    XR <- int.r$RR
  }
  if(multipl == "OR") {
    XR <- int.r$OR
  }

  ## reference level
  if(ncol(labels.A1) == 1) {
    ref.A1 <- labels.A1[rowSums(labels.A1) == 0, ncol(labels.A1)]
  }
  if(ncol(labels.A1) > 1) {
    ref.A1 <- labels.A1[rowSums(labels.A1[,-ncol(labels.A1)]) == 0, ncol(labels.A1)]
  }

  if(ncol(labels.A2) == 1) {
    ref.A2 <- labels.A2[rowSums(labels.A2) == 0, ncol(labels.A2)]
  }
  if(ncol(labels.A2) > 1) {
    ref.A2 <- labels.A2[rowSums(labels.A2[,-ncol(labels.A2)]) == 0, ncol(labels.A2)]
  }

  ##  table of marginal effects (with A2 or V nodes displayed by column)
  out.table <- data.frame(matrix("", nrow = (length(int.r$A1nodes) + 1) + (length(int.r$A1nodes) * 2),
                                 ncol = (length(A2nodes) + 1) + (length(A2nodes) * 2)) )
  names.A1 <- paste0(names(labels.A1)[ncol(labels.A1)],"=",labels.A1[,ncol(labels.A1)])
  if(ncol(labels.A1) == 1) {
    names.RD.A1 <- paste0("RD.",names(labels.A1)[ncol(labels.A1)],"=",
                          labels.A1[rowSums(labels.A1)!=0, ncol(labels.A1)]," vs.",ref.A1,
                          "|",names(labels.A2)[ncol(labels.A2)])
    names.RR.A1 <- paste0(multipl,".",names(labels.A1)[ncol(labels.A1)],"=",
                          labels.A1[rowSums(labels.A1)!=0,ncol(labels.A1)]," vs.",ref.A1,
                          "|",names(labels.A2)[ncol(labels.A2)])
  }
  if(ncol(labels.A1) > 1) {
    names.RD.A1 <- paste0("RD.",names(labels.A1)[ncol(labels.A1)],"=",
                          labels.A1[rowSums(labels.A1[,1:(ncol(labels.A1)-1)])!=0,
                                    ncol(labels.A1)]," vs.",ref.A1,
                          "|",names(labels.A2)[ncol(labels.A2)])
    names.RR.A1 <- paste0(multipl,".",names(labels.A1)[ncol(labels.A1)],"=",
                          labels.A1[rowSums(labels.A1[,1:(ncol(labels.A1)-1)])!=0,
                                    ncol(labels.A1)]," vs.",ref.A1,
                          "|",names(labels.A2)[ncol(labels.A2)])
  }

  names.A2 <- paste0(names(labels.A2)[ncol(labels.A2)],"=",labels.A2[,ncol(labels.A2)])
  if(ncol(labels.A2) == 1) {
    names.RD.A2 <- paste0("RD.",names(labels.A2)[ncol(labels.A2)],"=",
                          labels.A2[rowSums(labels.A2)!=0, ncol(labels.A2)]," vs.",ref.A2,
                          "|",names(labels.A1)[ncol(labels.A1)])
    names.RR.A2 <- paste0(multipl,".",names(labels.A2)[ncol(labels.A2)],"=",
                          labels.A2[rowSums(labels.A2)!=0, ncol(labels.A2)]," vs.",ref.A2,
                          "|",names(labels.A1)[ncol(labels.A1)])
  }

  if(ncol(labels.A2) > 1) {
    names.RD.A2 <- paste0("RD.",names(labels.A2)[ncol(labels.A2)],"=",
                          labels.A2[rowSums(labels.A2[,1:(ncol(labels.A2)-1)])!=0,
                                    ncol(labels.A2)]," vs.",ref.A2,
                          "|",names(labels.A1)[ncol(labels.A1)])
    names.RR.A2 <- paste0(multipl,".",names(labels.A2)[ncol(labels.A2)],"=",
                          labels.A2[rowSums(labels.A2[,1:(ncol(labels.A2)-1)])!=0,
                                    ncol(labels.A2)]," vs.",ref.A2,
                          "|",names(labels.A1)[ncol(labels.A1)])
  }

  names(out.table) <-c(names.A2, names.RD.A2, names.RR.A2)
  rownames(out.table) <- c(names.A1, names.RD.A1, names.RR.A1)

  int_res = int.r$int

  # p
  for(r in 1:nrow(labels.A1)) {
    for(c in 1:nrow(labels.A2)) {
      if(ncol(labels.A1) == 1) {
        index.A1 <- int.r$probs[,names(labels.A1)] == labels.A1[r,]
      }
      if(ncol(labels.A1) > 1) {
        index.A1 <- apply(apply(int.r$probs[,names(labels.A1)[-ncol(labels.A1)]], 1, function(row) row == labels.A1[r,-ncol(labels.A1)]),2,prod) == 1
      }
      if(ncol(labels.A2) == 1) {
        index.A2 <- int.r$probs[,names(labels.A2)] == labels.A2[c,]
      }
      if(ncol(labels.A2) > 1) {
        index.A2 <- apply(apply(int.r$probs[,names(labels.A2)[-ncol(labels.A2)]], 1, function(row) row == labels.A2[c,-ncol(labels.A2)]),2,prod) == 1
      }
      p <- int.r$probs[index.A1 & index.A2,"p"]
      p.lo <- int.r$probs[index.A1 & index.A2,"p.lo"]
      p.up <- int.r$probs[index.A1 & index.A2,"p.up"]

      out.table[r,c] <- paste0(round(p, digits = probab.digits),
                               " [",
                               round(p.lo, digits = probab.digits),
                               ",",
                               round(p.up, digits = probab.digits),
                               "]")
      rm(index.A1, index.A2, p, p.lo, p.up)
    }
  }

  # RD.A1
  for(i in 1:ncol(int.r$values.A1)) {
    index.col <- length(int.r$A1nodes) + length(A2nodes) + 4 * (i - 1) + 1
    for(j in 1:nrow(labels.A2)) {
      index.A1 <- int.r$RD[,colnames(int.r$values.A1)[i]] == 1
      if(ncol(labels.A2) == 1) {
        index.A2 <- int.r$RD[,names(labels.A2)] == labels.A2[j,]
      }
      if(ncol(labels.A2) > 1) {
        index.A2 <- apply(apply(int.r$RD[,names(labels.A2)[-ncol(labels.A2)]], 1, function(row) row == labels.A2[j,-ncol(labels.A2)]),2,prod) == 1
      }

      RD <- int.r$RD[index.A1 & index.A2, index.col]
      lo.RD <- int.r$RD[index.A1 & index.A2, index.col + 2]
      up.RD <- int.r$RD[index.A1 & index.A2, index.col + 3]

      out.table[nrow(int.r$values.A1) + i, j] <- paste0(round(RD, digits = probab.digits),
                               " [",
                               round(lo.RD, digits = probab.digits),
                               ",",
                               round(up.RD, digits = probab.digits),
                               "]")
      rm(index.A1, index.A2, RD, lo.RD, up.RD)
    }
  }

  # RD.A2
  for(i in 1:ncol(values.A2)) {
    index.col <- length(int.r$A1nodes) + length(A2nodes) + 4 * length(int.r$A1nodes) + 4 * (i - 1) + 1
    for(j in 1:nrow(labels.A1)) {
      if(ncol(labels.A1) == 1) {
        index.A1 <- int.r$RD[,names(labels.A1)] == labels.A1[j,]
      }
      if(ncol(labels.A1) > 1) {
        index.A1 <- apply(apply(int.r$RD[,names(labels.A1)[-ncol(labels.A1)]], 1, function(row) row == labels.A1[j,-ncol(labels.A1)]),2,prod) == 1
      }
      index.A2 <- int.r$RD[,colnames(values.A2)[i]] == 1

      RD <- int.r$RD[index.A1 & index.A2, index.col]
      lo.RD <- int.r$RD[index.A1 & index.A2, index.col + 2]
      up.RD <- int.r$RD[index.A1 & index.A2, index.col + 3]

      out.table[j, nrow(values.A2) + i] <- paste0(round(RD, digits = probab.digits),
                                                        " [",
                                                        round(lo.RD, digits = probab.digits),
                                                        ",",
                                                        round(up.RD, digits = probab.digits),
                                                        "]")
      rm(index.A1, index.A2, RD, lo.RD, up.RD)
    }
  }

  # XR.A1
  for(i in 1:ncol(int.r$values.A1)) {
    index.col <- length(int.r$A1nodes) + length(A2nodes) + 4 * (i - 1) + 1
    for(j in 1:nrow(labels.A2)) {
      index.A1 <- XR[,colnames(int.r$values.A1)[i]] == 1
      if(ncol(labels.A2) == 1) {
        index.A2 <- XR[,names(labels.A2)] == labels.A2[j,]
      }
      if(ncol(labels.A2) > 1) {
        index.A2 <- apply(apply(XR[,names(labels.A2)[-ncol(labels.A2)]], 1, function(row) row == labels.A2[j,-ncol(labels.A2)]),2,prod) == 1
      }

      XR.est <- XR[index.A1 & index.A2, index.col]
      lo.XR <- XR[index.A1 & index.A2, index.col + 2]
      up.XR <- XR[index.A1 & index.A2, index.col + 3]

      out.table[nrow(int.r$values.A1) + ncol(int.r$values.A1) + i, j] <- paste0(round(XR.est, digits = RR.OR.digits),
                                                                                " [",
                                                                                round(lo.XR, digits = RR.OR.digits),
                                                                                ",",
                                                                                round(up.XR, digits = RR.OR.digits),
                                                                                "]")
      rm(index.A1, index.A2, XR.est, lo.XR, up.XR)
    }
  }

  # XR.A2
  for(i in 1:ncol(values.A2)) {
    index.col <- length(int.r$A1nodes) + length(A2nodes) + 4 * length(int.r$A1nodes) + 4 * (i - 1) + 1
    for(j in 1:nrow(labels.A1)) {
      if(ncol(labels.A1) == 1) {
        index.A1 <- XR[,names(labels.A1)] == labels.A1[j,]
      }
      if(ncol(labels.A1) > 1) {
        index.A1 <- apply(apply(XR[,names(labels.A1)[-ncol(labels.A1)]], 1, function(row) row == labels.A1[j,-ncol(labels.A1)]),2,prod) == 1
      }
      index.A2 <- XR[,colnames(values.A2)[i]] == 1

      XR.est <- XR[index.A1 & index.A2, index.col]
      lo.XR <- XR[index.A1 & index.A2, index.col + 2]
      up.XR <- XR[index.A1 & index.A2, index.col + 3]

      out.table[j, nrow(values.A2) + ncol(values.A2) + i] <- paste0(round(XR.est, digits = RR.OR.digits),
                                                  " [",
                                                  round(lo.XR, digits = RR.OR.digits),
                                                  ",",
                                                  round(up.XR, digits = RR.OR.digits),
                                                  "]")
      rm(index.A1, index.A2, XR.est, lo.XR, up.XR)
    }
  }


  # interaction
  interaction.effects <- NULL
  for(i in 1:ncol(int.r$values.A1)) {
    for(j in 1:ncol(values.A2)) {
      index.A1 <- int_res[,colnames(int.r$values.A1)[i]] == 1
      index.A2 <- int_res[,colnames(values.A2)[j]] == 1
      lab.A1 <- labels.A1[labels.A1[,i] == 1,ncol(labels.A1)]
      lab.A2 <- labels.A2[labels.A2[,j] == 1,ncol(labels.A2)]

      interaction.effects <- c(interaction.effects,
                               paste0(colnames(labels.A1)[ncol(labels.A1)],"=",lab.A1," vs.",ref.A1," & ",colnames(labels.A2)[ncol(labels.A2)],"=",lab.A2," vs.",ref.A2),
                               paste0(". additive Interaction = ",
                                      round(int_res[index.A1 & index.A2,"a.INT"], digits = probab.digits),
                                      " [",
                                      round(int_res[index.A1 & index.A2,"lo.a.INT"], digits = probab.digits),
                                      ";",
                                      round(int_res[index.A1 & index.A2,"up.a.INT"], digits = probab.digits),
                                      "]"),
                               paste0(". RERI = ",
                                      round(int_res[index.A1 & index.A2,"RERI"], digits = RR.OR.digits),
                                      " [",
                                      round(int_res[index.A1 & index.A2,"lo.RERI"], digits = RR.OR.digits),
                                      ";",
                                      round(int_res[index.A1 & index.A2,"up.RERI"], digits = RR.OR.digits),
                                      "]"),
                               paste0(". multiplicative Interaction = ",
                                      round(int_res[index.A1 & index.A2,"m.INT"], digits = RR.OR.digits),
                                      " [",
                                      round(int_res[index.A1 & index.A2,"lo.m.INT"], digits = RR.OR.digits),
                                      ";",
                                      round(int_res[index.A1 & index.A2,"up.m.INT"], digits = RR.OR.digits),
                                      "]"))
    }
  }



  if (int.r$transformOutcome == TRUE) {
    ncol.without.RR <- nrow(labels.A2) + ncol(labels.A2[,-ncol(labels.A2)])
    nrow.without.RR <- nrow(labels.A1) + ncol(labels.A1[,-ncol(labels.A1)])
    out.table <- out.table[1:nrow.without.RR,1:ncol.without.RR]
  }

  return(list(out.table = out.table,
              interaction.effects = interaction.effects))
}




#' Create plot of the marginal interaction effects
#'
#' @param int.r an \code{int.r} data frame from the output obtained by the \code{estim.int.effects} function
#'
#' @return A plot of the marginal interaction effects using \code{ggplot2} package.
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
#' # Estimate quantities of interest for the interaction effect of A1 * A2 -> Y
#' interaction.det <- estim.int.effects(interaction.ltmle, estimator = "tmle")
#'
#' # Show results in a plot
#' plot <- out.int.fig(interaction.det)
out.int.fig <- function(int.r = int.r,
                        labels.A1 = NULL,
                        labels.A2 = NULL,
                        labels.V = NULL) {
  if(is.null(int.r$A2nodes)) {
    A2nodes <- int.r$Vnodes
    labels.A2 <- labels.V
    values.A2 <- int.r$values.V
  }
  if(is.null(int.r$Vnodes)) {
    A2nodes <- int.r$A2nodes
    values.A2 <- int.r$values.A2
  }

  probs <- merge(int.r$probs, labels.A1, sort = FALSE)
  probs <- merge(probs, labels.A2, sort = FALSE)

  RD <- merge(int.r$RD, labels.A1, sort = FALSE)
  RD <- merge(RD, labels.A2, sort = FALSE)

  RR <- merge(int.r$RR, labels.A1, sort = FALSE)
  RR <- merge(RR, labels.A2, sort = FALSE)

  OR <- merge(int.r$OR, labels.A1, sort = FALSE)
  OR <- merge(OR, labels.A2, sort = FALSE)

  probs[,ncol(probs)] <- as.factor(probs[,ncol(probs)])
  RD[,ncol(RD)] <- as.factor(RD[,ncol(RD)])
  RR[,ncol(RR)] <- as.factor(RR[,ncol(RR)])
  OR[,ncol(OR)] <- as.factor(OR[,ncol(OR)])

  probs[,ncol(probs) - 1] <- as.factor(probs[,ncol(probs) - 1])
  RD[,ncol(RD) - 1] <- as.factor(RD[,ncol(RD) - 1])
  RR[,ncol(RR) - 1] <- as.factor(RR[,ncol(RR) - 1])
  OR[,ncol(OR) - 1] <- as.factor(OR[,ncol(OR) - 1])

  for(c in 1:ncol(RD)) {
    RD[which(is.na(RD[,c])),c] <- 0
  }
  for(c in 1:ncol(RR)) {
    RR[which(is.na(RR[,c])),c] <- 1
  }
  for(c in 1:ncol(OR)) {
    OR[which(is.na(OR[,c])),c] <- 1
  }

  A1 <- A2 <- NULL
  RD.A1 <- RD.A1.lo <- RD.A1.up <- RD.A2 <- RD.A2.lo <- RD.A2.up <- NULL
  RR.A1 <- RR.A1.lo <- RR.A1.up <- RR.A2 <- RR.A2.lo <- RR.A2.up <- NULL
  p <- p.lo <- p.up <- NULL

  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- ggplot2::position_dodge(0.04) # move them .05 to the left and right

  g1 <- ggplot2::ggplot(probs) +
    ggplot2::aes(x = probs[,ncol(probs) - 1], color = probs[,ncol(probs)], y = p) + # A1, A2
    ggplot2::geom_line(ggplot2::aes(group = probs[,ncol(probs)]), linetype = 2, position = pd) + # A2
    ggplot2::geom_errorbar(ggplot2::aes(group = probs[,ncol(probs)], ymin = p.lo, ymax = p.up), # A2
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::xlab(names(probs)[ncol(probs) - 1]) + # A1
    ggplot2::ylab(int.r$Ynodes) +
    ggplot2::labs(title = paste0("Effect of '",names(probs)[ncol(probs) - 1], "' in groups of '",names(probs)[ncol(probs)],"'")) + # A1, A2
    ggplot2::guides(color = ggplot2::guide_legend(title = names(probs)[ncol(probs)]))

  g2 <- ggplot2::ggplot(probs) +
    ggplot2::aes(x = probs[,ncol(probs)], color = probs[,ncol(probs) - 1], y = p) + # A2, A1
    ggplot2::geom_line(ggplot2::aes(group = probs[,ncol(probs) - 1]), linetype = 2, position = pd) + # A1
    ggplot2::geom_errorbar(ggplot2::aes(group = probs[,ncol(probs) - 1], ymin = p.lo, ymax = p.up), # A1
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::xlab(names(probs)[ncol(probs)]) + # A2
    ggplot2::ylab(int.r$Ynodes) +
    ggplot2::labs(title = paste0("Effect of '",names(probs)[ncol(probs)], "' in groups of '",names(probs)[ncol(probs) - 1],"'")) + # A2, A1
    ggplot2::guides(color = ggplot2::guide_legend(title = names(probs)[ncol(probs) - 1]))

  # les suivants sont plus difficiles, il faudrait regrouper les RD.A1 dans une même colonne. Idem pour les RD.A2
  g3 <- ggplot2::ggplot(RD) +
    ggplot2::aes(x = A1, color = A2, y= RD.A1) + # A1, A2
    ggplot2::geom_line(ggplot2::aes(group = A2, linetype = A2)) + # A2
    ggplot2::geom_errorbar(ggplot2::aes(group = A2, ymin = RD.A1.lo, ymax = RD.A1.up), # A2, RD.A1.lo, RD.A1.up
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::xlab(int.r$Anodes[1]) +
    ggplot2::ylab(paste0("RD.",int.r$Anodes[1])) +
    ggplot2::labs(subtitle = "(RD scale)") +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[2])),
                    linetype = ggplot2::guide_legend(title = paste0(int.r$Anodes[2])))

  g4 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A2, color = A1, y= RD.A2) +
    ggplot2::geom_line(ggplot2::aes(group = A1, linetype = A1)) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A1, ymin = RD.A2.lo, ymax = RD.A2.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::xlab(int.r$Anodes[2]) +
    ggplot2::ylab(paste0("RD.",int.r$Anodes[2])) +
    ggplot2::labs(subtitle = "(RD scale)") +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[1])),
                    linetype = ggplot2::guide_legend(title = paste0(int.r$Anodes[1])))

  g5 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A1, color = A2, y = RR.A1) +
    ggplot2::geom_line(ggplot2::aes(group = A2, linetype = A2)) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A2, ymin = RR.A1.lo, ymax = RR.A1.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::xlab(int.r$Anodes[1]) +
    ggplot2::ylab(paste0("RR.",int.r$Anodes[1])) +
    ggplot2::labs(subtitle = "(RR scale)") +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[2])),
                    linetype = ggplot2::guide_legend(title = paste0(int.r$Anodes[2])))

  g6 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A2, color = A1, y = RR.A2) +
    ggplot2::geom_line(ggplot2::aes(group = A1, linetype = A1)) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A1, ymin = RR.A2.lo, ymax = RR.A2.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::xlab(int.r$Anodes[2]) +
    ggplot2::ylab(paste0("RR.",int.r$Anodes[2])) +
    ggplot2::labs(subtitle = "(RR scale)") +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[1])),
                    linetype = ggplot2::guide_legend(title = paste0(int.r$Anodes[1])))

  ggpubr::ggarrange(g1, g2, g3, g4, g5, g6, ncol = 2, nrow = 3)

  if (int.r$transformOutcome == TRUE) {
    ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow = 2)
  }
}
