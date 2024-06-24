#' Fitting a Marginal Structural Model to estimate marginal interaction effects
#'
#' \code{int.ltmleMSM} is used to fit a MSM using the ltmleMSM function, from the
#' \code{ltmle} package.
#'
#' Details to detail..
#'
#' @param data data frame following the time-ordering of the nodes. See help of
#' \code{ltmle} package.
#' @param A1nodes column names or indicies in \code{data} of treatment nodes for the 1st exposure
#' \code{c("A1.1","A1.2",...)}. Dummy variables for categorical exposures with > 2 levels.
#' @param A2nodes column names or indicies in \code{data} of treatment nodes for the 2nd exposure
#' \code{c("A2.1","A2.2",...)}. Dummy variables for categorical exposures with > 2 levels.
#' @param Vnodes column name or index in \code{data} of effect modifier node
#' \code{c("V.1","V.2",...)}. Dummy variables for categorical exposures with > 2 levels.
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
#' @return \code{int.ltmleMSM} returns a list 5 objects:
#'                 \itemize{ \item \code{ltmle_MSM} the output from the ltmleMSM function
#'                           \item \code{data} the data frame where the exposures names are \code{c(A1,A2)} and the outcome name is \code{Y}
#'                           \item \code{Anodes} the vector c("A1","A2") of exposure nodes
#'                           \item \code{Ynodes} the outcome variable Y
#'                           \item \code{bootstrap.res} a data frame containing estimation from bootstrap samples (for g-computation only)}
#' @export
#'
#' @examples
#' ## Example 1
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
#'# Estimate MSM parameters by g-computation
#' interaction.gcomp <- int.ltmleMSM(data = df,
#'                                   Qform = Q_formulas,
#'                                   gform = g_formulas,
#'                                   A1nodes = c("behav.2","behav.3"),
#'                                   A2nodes = c("env.2","env.3"),
#'                                   Vnodes = NULL,
#'                                   Lnodes = NULL,
#'                                   Ynodes = c("hlth.outcome"),
#'                                   SL.library = list(Q="SL.lm", g="SL.mean"),
#'                                   gcomp = TRUE,
#'                                   iptw.only = FALSE,
#'                                   survivalOutcome = FALSE,
#'                                   variance.method = "ic",
#'                                   B = 5, # it should be at least 1000 or 2000
#'                                   boot.seed = 54321)
#'
#' ## Example 2 - c(env.2, env3) are effect modifiers among baseline confounders
#' set.seed(12345)
#' df <- generate.data.multcat(N = 1000, b = param.causal.model.multcat())
#' head(df)
#' df <- data.frame(df[,c("conf1","conf2","conf3")],
#'                  env.2 = ifelse(df$env == 2, 1, 0),
#'                  env.3 = ifelse(df$env == 3, 1, 0),
#'                  behav.2 = ifelse(df$behav == 2, 1, 0),
#'                  behav.3 = ifelse(df$behav == 3, 1, 0),
#'                  hlth.outcome = df$hlth.outcome)
#' head(df)
#'
#' # Define Q and g formulas
#' # an (A1 * Effect.modifier) interaction term is recommended in the Q formula for the estimation
#' # of interaction effects
#' Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + behav.2 * env.2 + behav.2 * env.3 + behav.3 * env.2 + behav.3 * env.3")
#' g_formulas = c("behav.2 ~ conf1 + conf2",
#'                "behav.3 ~ conf1 + conf2 + behav.2") #c(behav2,behav.3) are the only exposure variables
#'
#' # Define SuperLearner libraries
#' SL.library = list(Q = list("SL.glm"), g = list("SL.glm"))
#'
#' # Estimate MSM parameters by IPTW and TMLE
#' interaction.ltmle <- int.ltmleMSM(data = df,
#'                                   Qform = Q_formulas,
#'                                   gform = g_formulas,
#'                                   A1nodes = c("behav.2","behav.3"),
#'                                   A2nodes = NULL,
#'                                   Vnodes = c("env.2","env.3"),
#'                                   Lnodes = NULL,
#'                                   Ynodes = c("hlth.outcome"),
#'                                   SL.library = SL.library,
#'                                   gcomp = FALSE,
#'                                   iptw.only = FALSE,
#'                                   survivalOutcome = FALSE,
#'                                   variance.method = "ic")
int.ltmleMSM <- function(data = data,
                         A1nodes = A1nodes, # c("A1.1", "A1.2",...)
                         A2nodes = A2nodes, # c("A2.1", "A2.2",...)
                         Vnodes = Vnodes,
                         Cnodes = NULL,
                         Lnodes = NULL, # c("L1", "L2", "L3")
                         Ynodes = Ynodes, # "Y"
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
                         gcomp = FALSE,
                         iptw.only = iptw.only,
                         deterministic.Q.function = NULL,
                         variance.method = "tmle",
                         observation.weights = NULL,
                         id = NULL,
                         B = 1000,
                         boot.seed = NULL) {

  # check arguments effect modifiers arguments
  if(!is.null(A2nodes)) {
    try(if(!is.null(Vnodes))
      stop("if A2nodes is not NULL, then Vnodes should be NULL"))
  }
  if(!is.null(Vnodes)) {
    try(if(!is.null(A2nodes))
      stop("if Vnode is not NULL, then A2nodes should be NULL"))
  }

  ## format data set
  baseline.Lnodes <- names(df)[1:(min(c(which(names(df) %in% A1nodes),
                                        which(names(df) %in% A2nodes))) - 1)]
  interm.Lnodes <- Lnodes[which(!(Lnodes %in% baseline.Lnodes))]

  if(!is.null(Vnodes)) {
    for(i in length(Vnodes)) {
      try(if(!(Vnodes[i] %in% baseline.Lnodes))
        stop(paste0("The effect modifier ",Vnodes[i], " should be in the set of baseline confounders")))
    }
  }

  # # NON EN FAIT IL VAUDRAIT MIEUX SAISIR DIRECTEMENT LES DUMMIES POUR EVITER UN TRAITEMENT DIFFERENT ENTRE EXPOSITION ET AUTRES VARIABLES +++
  # # AU MOMENT DE L'IMPORT DES DONNEES
  # # create dummy variables for the exposures (or effect modifier)
  # if(length(Anodes) == 2) {
  #   levels.A1 <- sort(unique(data[,Anodes[1]]))
  #   A1 <- list()
  #   index <- 1
  #   for(i in which(!(levels.A1 %in% Anodes.ref[1]))) {
  #     A1[[index]] <- ifelse(data[,Anodes[1]] == i, 1, 0)
  #     names(A1)[index] <- paste0(Anodes[1],".",i)
  #     index <- index + 1
  #   }
  #
  #   levels.A2 <- sort(unique(data[,Anodes[2]]))
  #   A2 <- list()
  #   index <- 1
  #   for(i in which(!(levels.A2 %in% Anodes.ref[2]))) {
  #     A2[[index]] <- ifelse(data[,Anodes[2]] == i, 1, 0)
  #     names(A2)[index] <- paste0(Anodes[2],".",i)
  #     index <- index + 1
  #   }
  #   rm(index)
  #
  #   data.int <- data[,baseline.Lnodes]
  #   for(i in 1:length(A1)) {
  #     data.int <- cbind(data.int, A1[[i]])
  #     names(data.int)[ncol(data.int)] <- names(A1)[i]
  #   }
  #   data.int <- cbind(data.int, data[,interm.Lnodes])
  #   for(i in 1:length(A1)) {
  #     data.int <- cbind(data.int, A2[[i]])
  #     names(data.int)[ncol(data.int)] <- names(A2)[i]
  #   }
  #   rm(i)
  #   data.int <- cbind(data.int, data[,Ynodes])
  #   names(data.int)[ncol(data.int)] <- Ynodes
  # }
  #
  #

  if(!is.null(A2nodes)) {
    Anodes <- c(A1nodes, A2nodes)
    nb.regimes <- (length(A1nodes) + 1) * (length(A2nodes) + 1)
    numAnodes <- length(A1nodes) + length(A2nodes)

    int.terms <- list()
    index <- 1
    for(i in 1:length(A1nodes)) {
      for(j in 1:length(A2nodes)) {
        int.terms[[index]] <- paste0(A1nodes[i]," * ",A2nodes[j])
        index <- index + 1
      }
    }

    working.msm <- paste0("Y ~")
    model.msm <- paste0(Ynodes, " ~")
    for(i in 1:length(int.terms)) {
      working.msm <- paste0(working.msm," + ",int.terms[[i]])
      model.msm <- paste0(model.msm," + ",int.terms[[i]])
    }
    rm(index,i,j,int.terms)

    # summary.measures = valeurs des coefficients du MSM associés à chaque régime
    # array: num.regimes x num.summary.measures x num.final.Ynodes -
    # measures summarizing the regimes that will be used on the right hand side of working.msm
    # (baseline covariates may also be used in the right hand side of working.msm and do not need to be included in summary.measures)
    summary.measures.reg <- array(NA, dim = c((length(A1nodes) + 1) * (length(A2nodes) + 1),
                                              (length(A1nodes) + length(A2nodes)),
                                              length(Ynodes)))
    # summary.measures.reg <- array(NA, dim = c((length(A1nodes) + 1) * (length(A2nodes) + 1),
    #                                           (length(A1nodes) + length(A2nodes) + length(A1nodes) * length(A2nodes)),
    #                                           length(Ynodes)))
    summary.measures.reg[,,1] <- unique(model.matrix(as.formula(model.msm), data = data))[,-1][,1:(length(A1nodes) + length(A2nodes))]
    colnames(summary.measures.reg) <- colnames(unique(model.matrix(as.formula(model.msm), data = data))[,-1][,1:(length(A1nodes) + length(A2nodes))])

    for(i in length(Anodes)) {
      try(if(!(Anodes[i] %in% colnames(summary.measures.reg)))
        stop("Problem in creation of the summary.measures argument"))
    }

    # regime=
    # binary array: n x numAnodes x numRegimes of counterfactual treatment or a list of 'rule' functions
    # colnames(regime) should be the same as Anodes (in the same order)
    regimes.MSM.temp <- array(NA, dim = c(nrow(data), numAnodes, nb.regimes))
    for(i in 1:nb.regimes) {
      for(r in 1:nrow(data)) {
        regimes.MSM.temp[r,,i] <- summary.measures.reg[i,1:numAnodes,1]
      }
    }
    colnames(regimes.MSM.temp) <- colnames(summary.measures.reg)[1:numAnodes]
    regimes.MSM <- array(NA, dim = c(nrow(data), numAnodes, nb.regimes))
    for(i in 1:length(Anodes)) {
      regimes.MSM[,i,] <- regimes.MSM.temp[,Anodes[i],]
    }
    colnames(regimes.MSM) <- Anodes
    rm(i,regimes.MSM.temp)

    try(if(dim(regimes.MSM)[3] != (length(A1nodes) + 1) * (length(A2nodes) + 1))
      stop("Interaction effects cannot be identified from this data set (MSM is not saturated)"))
  }

  if(!is.null(Vnodes)) {
    Anodes <- c(A1nodes)
    nb.regimes <- (length(A1nodes) + 1)
    numAnodes <- length(A1nodes)

    int.terms <- list()
    index <- 1
    for(i in 1:length(A1nodes)) {
      for(j in 1:length(Vnodes)) {
        int.terms[[index]] <- paste0(A1nodes[i]," * ",Vnodes[j])
        index <- index + 1
      }
    }

    working.msm <- paste0("Y ~")
    model.msm <- paste0(Ynodes, " ~")
    for(i in 1:length(int.terms)) {
      working.msm <- paste0(working.msm," + ",int.terms[[i]])
      model.msm <- paste0(model.msm," + ",int.terms[[i]])
    }
    rm(index,i,j,int.terms)

    # summary.measures = valeurs des coefficients du MSM associés à chaque régime
    # array: num.regimes x num.summary.measures x num.final.Ynodes -
    # measures summarizing the regimes that will be used on the right hand side of working.msm
    # (baseline covariates may also be used in the right hand side of working.msm and do not need to be included in summary.measures)
    summary.measures.reg <- array(0, dim = c(length(A1nodes) + 1,
                                              length(A1nodes),
                                              length(Ynodes)))
    for(i in 1:length(A1nodes)) {
      summary.measures.reg[1 + i,i,] <- 1
    }
    colnames(summary.measures.reg) <- A1nodes

    # regime=
    # binary array: n x numAnodes x numRegimes of counterfactual treatment or a list of 'rule' functions
    # colnames(regime) should be the same as Anodes (in the same order)
    regimes.MSM <- array(NA, dim = c(nrow(data), numAnodes, nb.regimes))
    for(i in 1:nb.regimes) {
      for(r in 1:nrow(data)) {
        regimes.MSM[r,,i] <- summary.measures.reg[i,1:numAnodes,1]
      }
    }
    colnames(regimes.MSM) <- A1nodes
  }

  # A1node.matrix <- matrix(0, nrow = (length(A1nodes) + 1), ncol = length(A1nodes))
  # for(i in 1:ncol(A1node.matrix)) {
  #   A1node.matrix[i+1,i] <- 1
  # }
  # colnames(A1node.matrix) <- A1nodes
  # A2node.matrix <- matrix(0, nrow = (length(A2nodes) + 1), ncol = length(A2nodes))
  # for(i in 1:ncol(A2node.matrix)) {
  #   A2node.matrix[i+1,i] <- 1
  # }
  # colnames(A2node.matrix) <- A2nodes
  #
  #
  #
  # Anode.matrix <- matrix(NA, nrow = 0, #nrow = nrow(A1node.matrix) * nrow(A2node.matrix),
  #                        ncol = length(A1nodes) * length(A2nodes))
  # for(i in 1:nrow(A1node.matrix)) {
  #   for(j in 1:nrow(A2node.matrix)) {
  #     Anode.matrix <- rbind(Anode.matrix,
  #                           c(A1node.matrix[i,],A2node.matrix[j,]))
  #   }
  # }
  # for(i in 1:length(A1nodes)) {
  #   for(j in (length(A1nodes) +1):(length(A1nodes) + length(A2nodes))) {
  #     Anode.matrix <- cbind(Anode.matrix,
  #                           as.data.frame(Anode.matrix)[,i] * as.data.frame(Anode.matrix)[,j])
  #   }
  # }







  if(gcomp == TRUE) {
    # test length SL.library$Q
    SL.library$Q <- ifelse(length(SL.library$Q) > 1, "glm", SL.library$Q) # +++++ DELETE +++++++++++ ?????

    # simplify SL.library$g because g functions are useless with g-computation
    SL.library$g <- "SL.mean"

    iptw.only <- FALSE
  }



  ltmle_MSM <- ltmle::ltmleMSM(data = data,
                               Anodes = Anodes,
                               Cnodes = Cnodes,
                               Lnodes = c(interm.Lnodes),
                               Ynodes = Ynodes,
                               survivalOutcome = survivalOutcome,
                               Qform = Qform,
                               gform = gform,
                               gbounds = gbounds,
                               Yrange = Yrange,
                               deterministic.g.function = deterministic.g.function,
                               SL.library = SL.library,
                               SL.cvControl = SL.cvControl,
                               regimes = regimes.MSM, # instead of abar
                               working.msm= working.msm,
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
  bootstrap.res <- NULL

  if(gcomp == TRUE) {
    # bootstrap.res <- data.frame("beta.Intercept" = rep(NA, B),
    #                             "beta.A1" = rep(NA, B),
    #                             "beta.A2" = rep(NA, B),
    #                             "beta.A1A2" = rep(NA, B))

    bootstrap.res.list <- list()

    if (is.null(boot.seed)) {
      print("boot.seed argument is null, please add a seed in the int.ltmleMSM function")}
    set.seed <- boot.seed

    for (b in 1:B){
      # sample the indices 1 to n with replacement
      bootIndices <- sample(1:nrow(data), replace=T)
      bootData <- data[bootIndices,]

      if ( round(b/100, 0) == b/100 ) print(paste0("bootstrap number ",b))

      suppressMessages(boot_ltmle_MSM <- ltmle::ltmleMSM(data = bootData,
                                                         Anodes = Anodes,
                                                         Cnodes = Cnodes,
                                                         Lnodes = c(interm.Lnodes),
                                                         Ynodes = Ynodes,
                                                         survivalOutcome = survivalOutcome,
                                                         Qform = Qform,
                                                         gform = gform,
                                                         gbounds = gbounds,
                                                         Yrange = Yrange,
                                                         deterministic.g.function = deterministic.g.function,
                                                         SL.library = SL.library,
                                                         SL.cvControl = SL.cvControl,
                                                         regimes = regimes.MSM, # instead of abar
                                                         working.msm= working.msm,
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
                                                         id = id))

      bootstrap.res.list[[b]] <- boot_ltmle_MSM$beta
    }

    bootstrap.res <- matrix(NA, nrow = B, ncol = length(bootstrap.res.list[[1]]))
    for(b in 1:B) {
      for(i in 1:length(bootstrap.res.list[[1]])) {
        bootstrap.res[b,i] <- bootstrap.res.list[[b]][i]
      }
    }
    colnames(bootstrap.res) <- names(boot_ltmle_MSM$beta)
    bootstrap.res <- data.frame(bootstrap.res)
  }

  return(list(ltmle_MSM = ltmle_MSM,
              working.msm = working.msm,
              data = data,
              A1nodes = A1nodes,
              A2nodes = A2nodes,
              Vnodes = Vnodes,
              Ynodes = Ynodes,
              bootstrap.res = bootstrap.res))
}
