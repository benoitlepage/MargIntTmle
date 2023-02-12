#' Fitting a Marginal Structural Model to estimate marginal interaction effects
#'
#' \code{int.ltmleMSM} is used to fit a MSM using the ltmleMSM function, from the \code{ltmle} package.
#'
#' Details to detail..
#'
#' @param data data frame following the time-ordering of the nodes. Should follow the format recommended for the \code{ltmle} package
#' @param Anodes column names in \code{data} for the two exposures \code{c(A1,A2)}
#' @param Cnodes used for censoring nodes in ltmleMSM function. \code{NULL} by default, survival is not yet implemented for the \code{int.ltmleMSM function}
#' @param Lnodes column names in \code{data} for confounders of the A1 -> Y and A2 -> Y relationships
#' @param Ynodes column names in \code{data} for the outcome node
#' @param survivalOutcome column names in \code{data} for the outcome node. \code{FALSE} by default, survival is not yet implemented for the \code{int.ltmleMSM function}
#' @param Qform
#' @param gform
#' @param gbounds
#' @param Yrange
#' @param deterministic.g.function
#' @param SL.library
#' @param SL.cvControl
#' @param final.Ynodes
#' @param stratify
#' @param msm.weights
#' @param estimate.time
#' @param gcomp
#' @param iptw.only
#' @param deterministic.Q.function
#' @param variance.method
#' @param observation.weights
#' @param id
#' @param B
#' @param boot.seed
#'
#' @return \code{int.ltmleMSM} returns 3 objects:
#'                 \itemize{ \item \code{ltmle_MSM} the output from the ltmleMSM function
#'                           \item \code{df.int} the data frame where the exposures names are \code{c(A1,A2)} and the outcome name is \code{Y}}
#'                           \item \code{bootstrap.res} a data frame containing estimation from bootstrap samples (for g-computation only)
#' @export
#'
#' @examples
#' set.seed(12345)
#' b = param.causal.model()
#' df <- generate.data(N = 10000, b = b)
#' summary(df)
#'
#' # Define Q and g formulas
#' # an A1 * A2 interaction term is recommended in the Q formula for the estimation of interaction effects
#' Q_formulas = c(Y="Q.kplus1 ~ conf1 + conf2 + conf3 + A1 * A2")
#' g_formulas = c("A1 ~ conf1 + conf2",
#'                "A2 ~ conf1 + conf3")
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
#'# Estimate MSM parameters by g-computation
#' interaction.gcomp <- int.ltmleMSM(data = df,
#'                                   Qform = Q_formulas,
#'                                   gform = g_formulas,
#'                                   Anodes = c("sex", "env"),
#'                                   Lnodes = c("conf1", "conf2", "conf3"),
#'                                   Ynodes = c("hlth.outcome"),
#'                                   SL.library = list(Q="SL.lm", g="SL.mean"),
#'                                   gcomp = TRUE,
#'                                   iptw.only = FALSE,
#'                                   survivalOutcome = FALSE,
#'                                   variance.method = "ic",
#'                                   B = 100,
#'                                   boot.seed = 54321)
int.ltmleMSM <- function(data = data,
                         Anodes = Anodes, # c(A1, A2)
                         Cnodes = NULL,
                         Lnodes = NULL, # c(L1, L2, L3)
                         Ynodes = Ynodes, # Y
                         survivalOutcome = FALSE,
                         Qform = Qform,
                         gform = gform,
                         gbounds = c(0.01, 1),
                         Yrange = NULL,
                         deterministic.g.function = NULL,
                         SL.library = list(Q="SL.glm",
                                           g="SL.glm"),
                         SL.cvControl = list(),
                         final.Ynodes = NULL, # Y
                         stratify = FALSE,
                         msm.weights = "empirical",
                         estimate.time = FALSE,
                         gcomp = gcomp,
                         iptw.only = iptw.only,
                         deterministic.Q.function = NULL,
                         variance.method = "ic",
                         observation.weights = NULL,
                         id = NULL,
                         B = 2000,
                         boot.seed = NULL) {

  require(ltmle)
  require(SuperLearner)

  df.int <- data.frame(data[,which(names(data) == Lnodes)],
                       A1 = data[,which(names(data) == Anodes[1])],
                       A2 = data[,which(names(data) == Anodes[2])],
                       Y = data[,which(names(data) == Ynodes)])

  # regime=
  # binary array: n x numAnodes x numRegimes of counterfactual treatment or a list of 'rule' functions
  regimes.MSM <- array(NA, dim = c(nrow(data), 2, 4)) # 2 variables d'exposition (A1, A2), 4 régimes d'exposition (0,0) (1,0) (0,1) (1,1)
  regimes.MSM[,,1] <- matrix(c(0,0), ncol = 2, nrow = nrow(data), byrow = TRUE) # exposé ni à A1, ni à A2
  regimes.MSM[,,2] <- matrix(c(1,0), ncol = 2, nrow = nrow(data), byrow = TRUE) # exposé à A1 uniquement
  regimes.MSM[,,3] <- matrix(c(0,1), ncol = 2, nrow = nrow(data), byrow = TRUE) # exposé à A2 uniquement
  regimes.MSM[,,4] <- matrix(c(1,1), ncol = 2, nrow = nrow(data), byrow = TRUE) # exposé à A1 et à A2

  # summary.measures = valeurs des coefficients du MSM associés à chaque régime
  # array: num.regimes x num.summary.measures x num.final.Ynodes -
  # measures summarizing the regimes that will be used on the right hand side of working.msm
  # (baseline covariates may also be used in the right hand side of working.msm and do not need to be included in summary.measures)
  summary.measures.reg <- array(NA, dim = c(4, 3, 1))
  summary.measures.reg[,,1] <- matrix(c(0, 0, 0, # aucun effet ni de A1, ni de A2
                                        1, 0, 0, # effet de A1 isolé
                                        0, 1, 0, # effet de A2 isolé
                                        1, 1, 1), # effet de A1 + A2 + A1:A2
                                      ncol = 3, nrow = 4, byrow = TRUE)
  colnames(summary.measures.reg) <- c("A1", "A2", "A1:A2")

  if(gcomp == TRUE) {
    # test length SL.library$Q
    SL.library$Q <- ifelse(length(SL.library$Q) > 1, "SL.glm", SL.library$Q)

    # simplify SL.library$g because g functions are useless with g-computation
    SL.library$g <- "SL.mean"

    iptw.only <- FALSE
  }

  ltmle_MSM <- ltmle::ltmleMSM(data = df.int,
                               Anodes = c("A1","A2"),
                               Cnodes = Cnodes,
                               Lnodes = Lnodes,
                               Ynodes = c("Y"),
                               survivalOutcome = survivalOutcome,
                               Qform = Qform,
                               gform = gform,
                               gbounds = gbounds,
                               Yrange = Yrange,
                               deterministic.g.function = deterministic.g.function,
                               SL.library = SL.library,
                               SL.cvControl = SL.cvControl,
                               regimes = regimes.MSM, # instead of abar
                               working.msm= "Y ~ A1 + A2 + A1:A2",
                               summary.measures = summary.measures.reg,
                               final.Ynodes = final.Ynodes,
                               stratify = stratify,
                               msm.weights = msm.weights,
                               estimate.time = estimate.time,
                               gcomp = gcomp,
                               iptw.only = iptw.only,
                               deterministic.Q.function = deterministic.Q.function,
                               variance.method = variance.method,
                               observation.weights = observation.weights,
                               id = id)

  bootstrap.res <- data.frame("beta.Intercept" = rep(NA, B),
                              "beta.A1" = rep(NA, B),
                              "beta.A2" = rep(NA, B),
                              "beta.A1A2" = rep(NA, B))

  if(gcomp == TRUE) {
    try(if(is.null(boot.seed))
      stop("boot.seed argument is null, please add a seed in the int.ltmleMSM function"))
    set.seed <- boot.seed

    for (b in 1:B){
      # sample the indices 1 to n with replacement
      bootIndices <- sample(1:nrow(data), replace=T)
      bootData <- df.int[bootIndices,]

      if ( round(b/100, 0) == b/100 ) print(paste0("bootstrap number ",b))

      boot_ltmle_MSM <- ltmleMSM(data = bootData,
                                 Anodes = c("A1","A2"),
                                 Lnodes = Lnodes,
                                 Ynodes = c("Y"),
                                 survivalOutcome = survivalOutcome,
                                 Qform = Q_formulas,
                                 gform = g_formulas,
                                 gbounds = gbounds,
                                 Yrange = Yrange,
                                 deterministic.g.function = deterministic.g.function,
                                 SL.library = SL.library,
                                 SL.cvControl = SL.cvControl,
                                 regimes = regimes.MSM, # instead of abar
                                 working.msm= "Y ~ A1 + A2 + A1:A2",
                                 summary.measures = summary.measures.reg,
                                 final.Ynodes = final.Ynodes,
                                 stratify = stratify,
                                 msm.weights = msm.weights,
                                 estimate.time = FALSE,
                                 gcomp = gcomp,
                                 iptw.only = iptw.only,
                                 deterministic.Q.function = deterministic.Q.function,
                                 variance.method = variance.method,
                                 observation.weights = observation.weights,
                                 id = id)

      bootstrap.res$beta.Intercept[b] <- boot_ltmle_MSM$beta["(Intercept)"]
      bootstrap.res$beta.A1[b] <- boot_ltmle_MSM$beta["A1"]
      bootstrap.res$beta.A2[b] <- boot_ltmle_MSM$beta["A2"]
      bootstrap.res$beta.A1A2[b] <- boot_ltmle_MSM$beta["A1:A2"]
    }
  }

  return(list(ltmle_MSM = ltmle_MSM,
              df.int = df.int,
              bootstrap.res = bootstrap.res))
}
