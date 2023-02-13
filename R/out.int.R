#' Create table of the marginal interaction effects
#'
#' @param int.r an output from \code{estim.int.effects} function
#'
#' @return
#' @export
#'
#' @examples
out.int.table <- function(int.r = int.r) {
  ##  tableau des effets marginaux
  out.table <- data.frame(c1 = rep("",4), c2 = rep("",4), c3 = rep("",4), c4 = rep("",4))
  names(out.table) <- c("A2=0", "A2=1", "RD.A2|A1", "RR.A2|A1")
  rownames(out.table) <- c("A1=0", "A1=1", "RD.A1|A2", "RR.A1|A2")
  # p
  out.table["A1=0","A2=0"] <- paste0("$p_{00}$=",round(int.r$p[which(int.r$A1==0 & int.r$A2==0)], digits = 3),
                                     " [",
                                     round(int.r$p.lo[which(int.r$A1==0 & int.r$A2==0)], digits = 3),
                                     ",",
                                     round(int.r$p.up[which(int.r$A1==0 & int.r$A2==0)], digits = 3),
                                     "]")
  out.table["A1=0","A2=1"] <- paste0("$p_{01}$=",round(int.r$p[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                     " [",
                                     round(int.r$p.lo[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                     ",",
                                     round(int.r$p.up[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                     "]")
  out.table["A1=1","A2=0"] <- paste0("$p_{10}$=",round(int.r$p[which(int.r$A1==1 & int.r$A2==0)], digits = 3),
                                     " [",
                                     round(int.r$p.lo[which(int.r$A1==1 & int.r$A2==0)], digits = 3),
                                     ",",
                                     round(int.r$p.up[which(int.r$A1==1 & int.r$A2==0)], digits = 3),
                                     "]")
  out.table["A1=1","A2=1"] <- paste0("$p_{11}$=",round(int.r$p[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                     " [",
                                     round(int.r$p.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                     ",",
                                     round(int.r$p.up[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                     "]")
  # RD
  out.table["A1=0","RD.A2|A1"] <- paste0(round(int.r$RD.A2[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                         " [",
                                         round(int.r$RD.A2.lo[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                         ",",
                                         round(int.r$RD.A2.up[which(int.r$A1==0 & int.r$A2==1)], digits = 3),
                                         "]")
  out.table["A1=1","RD.A2|A1"] <- paste0(round(int.r$RD.A2[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                         " [",
                                         round(int.r$RD.A2.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                         ",",
                                         round(int.r$RD.A2.up[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                         "]")
  out.table["RD.A1|A2","A2=0"] <- paste0(round(int.r$RD.A1[which(int.r$A1==1 & int.r$A2==0)],digits = 3),
                                         " [",
                                         round(int.r$RD.A1.lo[which(int.r$A1==1 & int.r$A2==0)],digits = 3),
                                         ",",
                                         round(int.r$RD.A1.up[which(int.r$A1==1 & int.r$A2==0)],digits = 3),
                                         "]")
  out.table["RD.A1|A2","A2=1"] <- paste0(round(int.r$RD.A1[which(int.r$A1==1 & int.r$A2==1)],digits = 3),
                                         " [",
                                         round(int.r$RD.A1.lo[which(int.r$A1==1 & int.r$A2==1)],digits = 3),
                                         ",",
                                         round(int.r$RD.A1.up[which(int.r$A1==1 & int.r$A2==1)],digits = 3),
                                         "]")

  # RR
  out.table["A1=0","RR.A2|A1"] <- paste0(round(int.r$RR.A2[which(int.r$A1==0 & int.r$A2==1)], digits = 2),
                                         " [",
                                         round(int.r$RR.A2.lo[which(int.r$A1==0 & int.r$A2==1)], digits = 2),
                                         ",",
                                         round(int.r$RR.A2.up[which(int.r$A1==0 & int.r$A2==1)], digits = 2),
                                         "]")
  out.table["A1=1","RR.A2|A1"] <- paste0(round(int.r$RR.A2[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                         " [",
                                         round(int.r$RR.A2.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                         ",",
                                         round(int.r$RR.A2.up[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                         "]")
  out.table["RR.A1|A2","A2=0"] <- paste0(round(int.r$RR.A1[which(int.r$A1==1 & int.r$A2==0)],digits = 2),
                                         " [",
                                         round(int.r$RR.A1.lo[which(int.r$A1==1 & int.r$A2==0)],digits = 2),
                                         ",",
                                         round(int.r$RR.A1.up[which(int.r$A1==1 & int.r$A2==0)],digits = 2),
                                         "]")
  out.table["RR.A1|A2","A2=1"] <- paste0(round(int.r$RR.A1[which(int.r$A1==1 & int.r$A2==1)],digits = 2),
                                         " [",
                                         round(int.r$RR.A1.lo[which(int.r$A1==1 & int.r$A2==1)],digits = 2),
                                         ",",
                                         round(int.r$RR.A1.up[which(int.r$A1==1 & int.r$A2==1)],digits = 2),
                                         "]")
  interaction.effects <- c(paste0("additive Interaction = ",
                                  round(int.r$a.INT[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                  " [",
                                  round(int.r$a.INT.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                  ";",
                                  round(int.r$a.INT.up[which(int.r$A1==1 & int.r$A2==1)], digits = 3),
                                  "]"),
                           paste0("RERI = ",
                                  round(int.r$RERI[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  " [",
                                  round(int.r$RERI.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  ";",
                                  round(int.r$RERI.up[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  "]"),
                           paste0("multiplicative Interaction = ",
                                  round(int.r$m.INT[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  " [",
                                  round(int.r$m.INT.lo[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  ";",
                                  round(int.r$m.INT.up[which(int.r$A1==1 & int.r$A2==1)], digits = 2),
                                  "]"))

  return(list(out.table = out.table,
              interaction.effects = interaction.effects))
}




#' Create plot of the marginal interaction effects
#'
#' @param int.r an output from \code{summary.int} function
#'
#' @return
#' @export
#'
#' @examples
out.int.fig <- function(int.r = int.r) {
  int.r$A1 <- as.factor(int.r$A1)
  int.r$A2 <- as.factor(int.r$A2)

  int.r$RD.A1[which(is.na(int.r$RD.A1))] <- 0
  int.r$RD.A2[which(is.na(int.r$RD.A2))] <- 0

  int.r$RR.A1[which(is.na(int.r$RR.A1))] <- 1
  int.r$RR.A2[which(is.na(int.r$RR.A2))] <- 1

  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- ggplot2::position_dodge(0.04) # move them .05 to the left and right


  g1 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A1, color = A2, y = p) +
    ggplot2::geom_line(aes(group = A2), linetype = 2, position = pd) +
    ggplot2::geom_errorbar(aes(group = A2, ymin = int.r$p.lo, ymax = int.r$p.up),
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::labs(title = "Effect of A1 in groups of A2")

  g2 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A2, color = A1, y = p) +
    ggplot2::geom_line(aes(group = A1), linetype = 2, position = pd) +
    ggplot2::geom_errorbar(aes(group = A1, ymin = int.r$p.lo, ymax = int.r$p.up),
                           width=.1, position = pd) +
    ggplot2::geom_point(position = pd) +
    ggplot2::labs(title = "Effect of A2 in groups of A1")

  g3 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A1, color = A2, y= RD.A1) +
    ggplot2::geom_line(aes(group = A2, linetype = A2)) +
    ggplot2::geom_errorbar(aes(group = A2, ymin = int.r$RD.A1.lo, ymax = int.r$RD.A1.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::labs(subtitle = "(RD scale)")

  g4 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A2, color = A1, y= RD.A2) +
    ggplot2::geom_line(aes(group = A1, linetype = A1)) +
    ggplot2::geom_errorbar(aes(group = A1, ymin = int.r$RD.A2.lo, ymax = int.r$RD.A2.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::labs(subtitle = "(RD scale)")

  g5 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A1, color = A2, y = RR.A1) +
    ggplot2::geom_line(aes(group = A2, linetype = A2)) +
    ggplot2::geom_errorbar(aes(group = A2, ymin = int.r$RR.A1.lo, ymax = int.r$RR.A1.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::labs(subtitle = "(RR scale)")

  g6 <- ggplot2::ggplot(int.r) +
    ggplot2::aes(x = A2, color = A1, y = RR.A2) +
    ggplot2::geom_line(aes(group = A1, linetype = A1)) +
    ggplot2::geom_errorbar(aes(group = A1, ymin = int.r$RR.A2.lo, ymax = int.r$RR.A2.up),
                           width=.1, position = pd) +
    ggplot2::geom_point() +
    ggplot2::labs(subtitle = "(RR scale)")

  g1 + g2 + g3 + g4 + g5 + g6 + plot_layout(ncol = 2, nrow = 3)
}
