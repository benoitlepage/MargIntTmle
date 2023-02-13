#' Fitting a Marginal Structural Model to estimate marginal interaction effects
#'
#' \code{int.ltmleMSM} is used to fit a MSM using the ltmleMSM function, from the
#' \code{ltmle} package.
#'
#' Details to detail..
#'
#' @param data data frame following the time-ordering of the nodes. See help of
#' \code{ltmle} package.
#' @param Anodes column names or indicies in \code{data} of treatment nodes
#' \code{c(A1,A2)}
#' @param Cnodes olumn names or indicies in \code{data} of censoring nodes.
#' \code{NULL} by default, survival is not yet implemented for the
#' \code{int.ltmleMSM function}
#' @param Lnodes column names or indicies in \code{data} of time-dependent
#' covariate nodes (confounders of the A1 -> Y and A2 -> Y
#' relationships)
#' @param Ynodes column names or indicies in \code{data} of outcome nodes
#' @param survivalOutcome If \code{TRUE}, then Y nodes are indicators of an
#' event, and if Y at some time point is 1, then all following should be 1.
#' Required to be \code{TRUE} or \code{FALSE} if outcomes are binary and there
#' are multiple Ynodes. \code{FALSE} by default, survival is not yet implemented
#' for the \code{int.ltmleMSM function}
#' @param Qform character vector of regression formulas for \eqn{\bar{Q}} function.
#' See 'Examples' and help of \code{ltmle} package.
#' @param gform character vector of regression formulas for gg or a matrix/array of
#' prob(A1=1) and prob(A2=1). See 'Examples' and help of \code{ltmle} package.
#' @param gbounds lower and upper bounds on estimated cumulative probabilities for
#' g-factors. Vector of length 2, order unimportant.
#' @param Yrange specify the range of all Y nodes. See 'Details'. See help of
#' \code{ltmle} package.
#' @param deterministic.g.function optional information on A and C nodes that are
#' given deterministically. See help of \code{ltmle} package. Default NULL
#' indicates no deterministic links. (? does not work with MSM ?)
#' @param SL.library optional character vector of libraries to pass to
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. \code{NULL} indicates
#' \link{glm} should be called instead of
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. '\code{default}'
#' indicates a standard set of libraries. May be separately specified for
#' \eqn{Q} and \eqn{g}. See help of \code{ltmle} package.
#' @param SL.cvControl optional list to be passed as \code{cvControl} to
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}
#' @param final.Ynodes vector subset of Ynodes - used in MSM to pool over a set
#' of outcome nodes.
#' @param stratify if TRUE stratify on following abar when estimating Q and g.
#' If FALSE, pool over abar.
#' @param msm.weights projection weights for the working MSM. If "empirical",
#' weight by empirical proportions of rows matching each regime for each
#' final.Ynode, with duplicate regimes given zero weight. If \code{NULL}, no
#' weights. Or an array of user-supplied weights with dimensions c(n,
#' num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes).
#' @param estimate.time if \code{TRUE}, run an initial estimate using only 50
#' observations and use this to print a very rough estimate of the total time
#' to completion. No action if there are fewer than 50 observations. \code{FALSE}
#' by default.
#' @param gcomp if \code{TRUE}, run the maximum likelihood based G-computation
#' estimate \emph{instead} of TMLE. 95% confidence intervals will be estimated
#' by boostratp
#' @param iptw.only by default (\code{iptw.only = FALSE}), both TMLE and IPTW
#' are run in \code{ltmleMSM}. If \code{iptw.only = TRUE},
#' only IPTW is run, which is faster.
#' @param deterministic.Q.function deterministic.Q.function optional information
#' on Q given deterministically. See help of \code{ltmle} package. Default
#' \code{NULL} indicates no deterministic links.
#' @param variance.method Method for estimating variance of TMLE.
#' One of "ic", "tmle", "iptw". If "tmle", compute both the robust variance
#' estimate using TMLE and the influence curve based variance estimate (use the
#' larger of the two). If "iptw", compute both the robust variance
#' estimate using IPTW and the influence curve based variance estimate (use the
#' larger of the two). If "ic", only compute the influence curve based
#' variance estimate. "ic" is fastest, but may be substantially
#' anti-conservative if there are positivity violations or rare outcomes. "tmle" is
#' slowest but most robust if there are positivity violations or rare outcomes.
#' "iptw" is a compromise between speed and robustness.
#' variance.method="tmle" or "iptw" are not yet available with non-binary outcomes,
#' gcomp=TRUE, stratify=TRUE, or deterministic.Q.function.
#' variance.method="tmle" or "iptw" are not available with gcomp=TRUE (only a
#' bootstrap method will be applied).
#' @param observation.weights observation (sampling) weights. Vector of length
#' n. If \code{NULL}, assumed to be all 1.
#' @param id Household or subject identifiers. Vector of length n or \code{NULL}.
#' Integer, factor, or character recommended, but any type that can be coerced
#' to factor will work. \code{NULL} means all distinct ids.
#' @param B if g-comp=TRUE, number of boostrap sample
#' @param boot.seed if g-comp=TRUE, seed for sampling bootstrap data sets
#'
#' @return \code{int.ltmleMSM} returns 3 objects:
#'                 \itemize{ \item \code{ltmle_MSM} the output from the ltmleMSM function
#'                           \item \code{df.int} the data frame where the exposures names are \code{c(A1,A2)} and the outcome name is \code{Y}
#'                           \item \code{bootstrap.res} a data frame containing estimation from bootstrap samples (for g-computation only)}
#' @export
#'
#' @examples
#' set.seed(12345)
#' b = param.causal.model()
#' df <- generate.data(N = 10000, b = b)
#' summary(df)
#'
#' # Define Q and g formulas
#' # an A1 * A2 interaction term is recommended in the Q formula for the estimation
#' # of interaction effects
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
#'                                   B = 5,
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
  # format data set with baseline confounders, 2 exposures (A1,A2) and the outcome Y
  # TO DO: modify to enable intermediate confounders between A1 and A2
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

      boot_ltmle_MSM <- ltmle::ltmleMSM(data = bootData,
                                        Anodes = c("A1","A2"),
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
