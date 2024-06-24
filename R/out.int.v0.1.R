## FIRST VERSION

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
                          RR.digits = 2) {
  ##  table of marginal effects
  out.table <- data.frame(c1 = rep("",4), c2 = rep("",4), c3 = rep("",4), c4 = rep("",4))
  names(out.table) <- c("A2=0", "A2=1", "RD.A2|A1", "RR.A2|A1")
  rownames(out.table) <- c("A1=0", "A1=1", "RD.A1|A2", "RR.A1|A2")

  int_res = int.r$int.r

  # p
  out.table["A1=0","A2=0"] <- paste0("$p_{00}$=",round(int_res$p[which(int_res$A1==0 & int_res$A2==0)], digits = probab.digits),
                                     " [",
                                     round(int_res$p.lo[which(int_res$A1==0 & int_res$A2==0)], digits = probab.digits),
                                     ",",
                                     round(int_res$p.up[which(int_res$A1==0 & int_res$A2==0)], digits = probab.digits),
                                     "]")
  out.table["A1=0","A2=1"] <- paste0("$p_{01}$=",round(int_res$p[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                     " [",
                                     round(int_res$p.lo[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                     ",",
                                     round(int_res$p.up[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                     "]")
  out.table["A1=1","A2=0"] <- paste0("$p_{10}$=",round(int_res$p[which(int_res$A1==1 & int_res$A2==0)], digits = probab.digits),
                                     " [",
                                     round(int_res$p.lo[which(int_res$A1==1 & int_res$A2==0)], digits = probab.digits),
                                     ",",
                                     round(int_res$p.up[which(int_res$A1==1 & int_res$A2==0)], digits = probab.digits),
                                     "]")
  out.table["A1=1","A2=1"] <- paste0("$p_{11}$=",round(int_res$p[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                     " [",
                                     round(int_res$p.lo[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                     ",",
                                     round(int_res$p.up[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                     "]")
  # RD
  out.table["A1=0","RD.A2|A1"] <- paste0(round(int_res$RD.A2[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                         " [",
                                         round(int_res$RD.A2.lo[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                         ",",
                                         round(int_res$RD.A2.up[which(int_res$A1==0 & int_res$A2==1)], digits = probab.digits),
                                         "]")
  out.table["A1=1","RD.A2|A1"] <- paste0(round(int_res$RD.A2[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                         " [",
                                         round(int_res$RD.A2.lo[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                         ",",
                                         round(int_res$RD.A2.up[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                         "]")
  out.table["RD.A1|A2","A2=0"] <- paste0(round(int_res$RD.A1[which(int_res$A1==1 & int_res$A2==0)],digits = probab.digits),
                                         " [",
                                         round(int_res$RD.A1.lo[which(int_res$A1==1 & int_res$A2==0)],digits = probab.digits),
                                         ",",
                                         round(int_res$RD.A1.up[which(int_res$A1==1 & int_res$A2==0)],digits = probab.digits),
                                         "]")
  out.table["RD.A1|A2","A2=1"] <- paste0(round(int_res$RD.A1[which(int_res$A1==1 & int_res$A2==1)],digits = probab.digits),
                                         " [",
                                         round(int_res$RD.A1.lo[which(int_res$A1==1 & int_res$A2==1)],digits = probab.digits),
                                         ",",
                                         round(int_res$RD.A1.up[which(int_res$A1==1 & int_res$A2==1)],digits = probab.digits),
                                         "]")

  # RR
  out.table["A1=0","RR.A2|A1"] <- paste0(round(int_res$RR.A2[which(int_res$A1==0 & int_res$A2==1)], digits = RR.digits),
                                         " [",
                                         round(int_res$RR.A2.lo[which(int_res$A1==0 & int_res$A2==1)], digits = RR.digits),
                                         ",",
                                         round(int_res$RR.A2.up[which(int_res$A1==0 & int_res$A2==1)], digits = RR.digits),
                                         "]")
  out.table["A1=1","RR.A2|A1"] <- paste0(round(int_res$RR.A2[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                         " [",
                                         round(int_res$RR.A2.lo[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                         ",",
                                         round(int_res$RR.A2.up[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                         "]")
  out.table["RR.A1|A2","A2=0"] <- paste0(round(int_res$RR.A1[which(int_res$A1==1 & int_res$A2==0)],digits = RR.digits),
                                         " [",
                                         round(int_res$RR.A1.lo[which(int_res$A1==1 & int_res$A2==0)],digits = RR.digits),
                                         ",",
                                         round(int_res$RR.A1.up[which(int_res$A1==1 & int_res$A2==0)],digits = RR.digits),
                                         "]")
  out.table["RR.A1|A2","A2=1"] <- paste0(round(int_res$RR.A1[which(int_res$A1==1 & int_res$A2==1)],digits = RR.digits),
                                         " [",
                                         round(int_res$RR.A1.lo[which(int_res$A1==1 & int_res$A2==1)],digits = RR.digits),
                                         ",",
                                         round(int_res$RR.A1.up[which(int_res$A1==1 & int_res$A2==1)],digits = RR.digits),
                                         "]")
  interaction.effects <- c(paste0("additive Interaction = ",
                                  round(int_res$a.INT[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                  " [",
                                  round(int_res$a.INT.lo[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                  ";",
                                  round(int_res$a.INT.up[which(int_res$A1==1 & int_res$A2==1)], digits = probab.digits),
                                  "]"),
                           paste0("RERI = ",
                                  round(int_res$RERI[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  " [",
                                  round(int_res$RERI.lo[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  ";",
                                  round(int_res$RERI.up[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  "]"),
                           paste0("multiplicative Interaction = ",
                                  round(int_res$m.INT[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  " [",
                                  round(int_res$m.INT.lo[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  ";",
                                  round(int_res$m.INT.up[which(int_res$A1==1 & int_res$A2==1)], digits = RR.digits),
                                  "]"))

  names(out.table) <- c(paste0(int.r$Anodes[2],"=0"),
                        paste0(int.r$Anodes[2],"=1"),
                        paste0("RD.",int.r$Anodes[2],"|",int.r$Anodes[1]),
                        paste0("RR.",int.r$Anodes[2],"|",int.r$Anodes[1]))
  rownames(out.table) <- c(paste0(int.r$Anodes[1],"=0"),
                           paste0(int.r$Anodes[1],"=1"),
                           paste0("RD.",int.r$Anodes[1],"|",int.r$Anodes[2]),
                           paste0("RR.",int.r$Anodes[1],"|",int.r$Anodes[2]))

  if (int.r$transformOutcome == TRUE) {
    out.table <- out.table[1:3,1:3]
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
out.int.fig <- function(int.r = int.r) {
  int_res = int.r$int.r

  int_res$A1 <- as.factor(int_res$A1)
  int_res$A2 <- as.factor(int_res$A2)

  int_res$RD.A1[which(is.na(int_res$RD.A1))] <- 0
  int_res$RD.A2[which(is.na(int_res$RD.A2))] <- 0

  int_res$RR.A1[which(is.na(int_res$RR.A1))] <- 1
  int_res$RR.A2[which(is.na(int_res$RR.A2))] <- 1

  A1 <- A2 <- NULL
  RD.A1 <- RD.A1.lo <- RD.A1.up <- RD.A2 <- RD.A2.lo <- RD.A2.up <- NULL
  RR.A1 <- RR.A1.lo <- RR.A1.up <- RR.A2 <- RR.A2.lo <- RR.A2.up <- NULL
  p <- p.lo <- p.up <- NULL

  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- ggplot2::position_dodge(0.04) # move them .05 to the left and right

  g1 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A1, color = A2, y = p) +
    ggplot2::geom_line(ggplot2::aes(group = A2), linetype = 2, position = pd) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A2, ymin = p.lo, ymax = p.up),
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::xlab(int.r$Anodes[1]) +
    ggplot2::ylab(int.r$Ynodes) +
    ggplot2::labs(title = paste0("Effect of '",int.r$Anodes[1], "' in groups of '",int.r$Anodes[2],"'")) +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[2])))

  g2 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A2, color = A1, y = p) +
    ggplot2::geom_line(ggplot2::aes(group = A1), linetype = 2, position = pd) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A1, ymin = p.lo, ymax = p.up),
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::xlab(int.r$Anodes[2]) +
    ggplot2::ylab(int.r$Ynodes) +
    ggplot2::labs(title = paste0("Effect of '",int.r$Anodes[2], "' in groups of '",int.r$Anodes[1],"'")) +
    ggplot2::guides(color = ggplot2::guide_legend(title = paste0(int.r$Anodes[1])))

  g3 <- ggplot2::ggplot(int_res) +
    ggplot2::aes(x = A1, color = A2, y= RD.A1) +
    ggplot2::geom_line(ggplot2::aes(group = A2, linetype = A2)) +
    ggplot2::geom_errorbar(ggplot2::aes(group = A2, ymin = RD.A1.lo, ymax = RD.A1.up),
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
