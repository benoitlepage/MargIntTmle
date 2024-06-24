#' Define parameters to simulate data
#'
#' This function defines parameters that are used to simulate illustrative data sets
#' consisting of three baseline confounders \eqn{L_1, L_2, L_3}, two exposures of interest
#' \eqn{A_1, A_2} and one outcome \eqn{Y}. All the variables are binary.
#'
#' @param p_L1 Probability \eqn{P(L_1 = 1) = p_{L_1}}.
#' @param p_L2 Probability \eqn{P(L_2 = 1) = p_{L_2}}.
#' @param p_L3 Probability \eqn{P(L_3 = 1) = p_{L_3}}.
#' @param b_A1 Probability \eqn{P(A_1 = 1|L_1 = 0, L_2 = 0) = \beta_{A_1}}.
#' @param b_L1_A1 Additive effect \eqn{\beta_{L_1,A_1}} of \eqn{L_1 \rightarrow A_1}.
#' @param b_L2_A1 Additive effect \eqn{\beta_{L_2,A_1}} of \eqn{L_2 \rightarrow A_1}.
#' @param b_A2 Probability \eqn{P(A_2 = 1|L_1 = 0, L_3 = 0) = \beta_{A_2}}.
#' @param b_L1_A2 Additive effect \eqn{\beta_{L_1,A_2}} of \eqn{L_1 \rightarrow A_2}.
#' @param b_L3_A2 Additive effect \eqn{\beta_{L_3,A_2}} of \eqn{L_3 \rightarrow A_2}.
#' @param Y_type Is Y "binary" or "continuous"?
#' @param b_Y Probability \eqn{P(Y = 1|L_1 = 0, L_2 = 0, L_3 = 0, A_1 = 0, A_2 = 0) = \beta_{Y}} if Y is binary;
#' or expectation \eqn{\mathbb{E}(Y|L_1 = 0, L_2 = 0, L_3 = 0, A_1 = 0, A_2 = 0) = \beta_{Y}}
#' @param b_L1_Y Additive effect \eqn{\beta_{L_1,Y}} of \eqn{L_1 \rightarrow Y}.
#' @param b_L2_Y Additive effect \eqn{\beta_{L_2,Y}} of \eqn{L_2 \rightarrow Y}.
#' @param b_L3_Y Additive effect \eqn{\beta_{L_3,Y}} of \eqn{L_3 \rightarrow Y}.
#' @param b_A1_Y Additive effect \eqn{\beta_{A_1,Y}} of \eqn{A_1 \rightarrow Y}.
#' @param b_A2_Y Additive effect \eqn{\beta_{A_2,Y}} of \eqn{A_2 \rightarrow Y}.
#' @param b_A1A2_Y Interaction additive effect \eqn{\beta_{A_1 \ast A_2,Y}} of \eqn{(A_1 \ast A_2) \rightarrow Y}.
#' @param se_Y Standard error of the outcome for quantitative outcomes, should be null for binary outcomes.
#'
#' @return A list of 4 vectors of parameters.
#' @export
#'
#' @examples
#' param.causal.model()
param.causal.model <- function(p_L1 = 0.50, p_L2 = 0.20, p_L3 = 0.70,  # baseline confounders
                               b_A1 = 0.10, b_L1_A1 = 0.15, b_L2_A1 = 0.25,  # exposure A1
                               b_A2 = 0.15, b_L1_A2 = 0.20, b_L3_A2 = 0.20,  # exposure A2
                               Y_type = "binary", # or "continuous"
                               b_Y = 0.10,   # outcome Y
                               b_L1_Y = 0.02,
                               b_L2_Y = 0.02,
                               b_L3_Y = -0.02,
                               b_A1_Y = 0.3,
                               b_A2_Y = 0.1,
                               b_A1A2_Y = 0.4,
                               se_Y = NULL) {  # b_A1A2_Y is the interaction effect A1 * A2 -> Y

  # check the sum of A1 parameters is not greater than 100%
  try(if(b_A1 + b_L1_A1 + b_L1_A1 > 1)
    stop("The sum of parameters to simulate A1 should not be greater than 100%"))

  # check the sum of A2 parameters is not greater than 100%
  try(if(b_A2 + b_L1_A2 + b_L3_A2 > 1)
    stop("The sum of parameters to simulate A2 should not be greater than 100%"))

  if (Y_type == "binary") {
    # check the sum of Y parameters is not greater than 100% nor less than 0%
    try(if(b_Y + b_L1_Y + b_L2_Y + b_L3_Y + b_A1_Y + b_A2_Y + b_A1A2_Y > 1)
      stop("The sum of parameters to simulate Y should not be greater than 100%"))
    try(if(b_Y + b_L1_Y + b_L2_Y + b_L3_Y + b_A1_Y + b_A2_Y + b_A1A2_Y < 0)
      stop("The sum of parameters to simulate Y should not be less than 0%"))
    try(if(!is.null(se_Y))
      stop("se_Y should be NULL for binary outcomes Y"))
  }
  if (Y_type == "continuous") {
    try(if(is.null(se_Y))
      stop("se_Y should not be NULL for continuous outcomes Y"))
  }

  coef <- list(c(p_L1 = p_L1, p_L2 = p_L2, p_L3 = p_L3),
               c(b_A1 = b_A1, b_L1_A1 = b_L1_A1, b_L2_A1 = b_L2_A1),
               c(b_A2 = b_A2, b_L1_A2 = b_L1_A2, b_L3_A2 = b_L3_A2),
               c(b_Y = b_Y, b_L1_Y = b_L1_Y, b_L2_Y = b_L2_Y, b_L3_Y = b_L3_Y,
                 b_A1_Y = b_A1_Y, b_A2_Y = b_A2_Y, b_A1A2_Y = b_A1A2_Y),
               c(se_Y = se_Y))
  return(coef)
}


#' Simulate an illustrative data set
#'
#' Simulate an illustrative data set of size \code{N} using the parameters \code{b}
#' defined using the \code{param.causal.model()} function.
#'
#' @param N Sample size of the illustrative data set.
#' @param b A list of 4 vectors of parameters,
#' @param Y_type Is Y "binary" or "continuous"?
#'
#'
#' @return A data frame of size \code{N} using the parameters \code{b},
#' generated from the following structural causal model:
#'
#' \eqn{P(L_1 = 1) = p_{L_1}}
#'
#' \eqn{P(L_2 = 1) = p_{L_2}}
#'
#' \eqn{P(L_3 = 1) = p_{L_3}}
#'
#' \eqn{P(A_1 = 1) = \beta_{A_1} + \beta_{L_1,A_1} L_1 + \beta_{L_2,A_1} L_2}
#'
#' \eqn{P(A_2 = 1) = \beta_{A_2} + \beta_{L_1,A_2} L_1 + \beta_{L_3,A_2} L_3}
#'
#' If Y is binary,
#' \eqn{P(Y = 1) = \beta_{Y} + \beta_{L_1,Y} L_1 + \beta_{L_2,Y} L_2 + \beta_{L_3,Y} L_3 + \beta_{A_1,Y} A_1 + \beta_{A_2,Y} A_2 + \beta_{A_1 \ast A_2,Y} (A_1 \ast A_2)}
#'
#' If Y is continuous,
#' \eqn{\mathbb{E}(Y|A_1,A_2,L) = \beta_{Y} + \beta_{L_1,Y} L_1 + \beta_{L_2,Y} L_2 + \beta_{L_3,Y} L_3 + \beta_{A_1,Y} A_1 + \beta_{A_2,Y} A_2 + \beta_{A_1 \ast A_2,Y} (A_1 \ast A_2)}
#'
#' The names of the baseline confounders \eqn{L_1, L_2, L_3} are \code{conf1}, \code{conf2}
#' and \code{conf3}. The names of the exposure \eqn{A_1} is \code{sex}, the name of
#' the exposure \eqn{A_2} is \code{env}. The name of the outcome \eqn{Y} is
#' \code{hlth.outcome}.
#'
#' @export
#'
#' @examples
#' set.seed(12345)
#' df <- generate.data(N = 1000, b = param.causal.model())
#' head(df)
generate.data <- function(N, b =  param.causal.model(), Y_type = "binary") {
  conf1 <- rbinom(N, size = 1, prob = b[[1]]["p_L1"])
  conf2 <- rbinom(N, size = 1, prob = b[[1]]["p_L2"])
  conf3 <- rbinom(N, size = 1, prob = b[[1]]["p_L3"])

  sex <- rbinom(N, size = 1, prob = b[[2]]["b_A1"] +
                  (b[[2]]["b_L1_A1"] * conf1) + (b[[2]]["b_L2_A1"] * conf2))
  env <- rbinom(N, size = 1, prob = b[[3]]["b_A2"] +
                  (b[[3]]["b_L1_A2"] * conf1) + (b[[3]]["b_L3_A2"] * conf3))

  if (Y_type == "binary") {
    hlth.outcome <- stats::rbinom(N, size = 1, prob = (b[[4]]["b_Y"] +
                                                         (b[[4]]["b_L1_Y"] * conf1) +
                                                         (b[[4]]["b_L2_Y"] * conf2) +
                                                         (b[[4]]["b_L3_Y"] * conf3) +
                                                         (b[[4]]["b_A1_Y"] * sex) +
                                                         (b[[4]]["b_A2_Y"] * env) +
                                                         (b[[4]]["b_A1A2_Y"] * sex * env)) )
  }
  if (Y_type == "continuous") {
    hlth.outcome <- stats::rnorm(N, mean = (b[[4]]["b_Y"] +
                                              (b[[4]]["b_L1_Y"] * conf1) +
                                              (b[[4]]["b_L2_Y"] * conf2) +
                                              (b[[4]]["b_L3_Y"] * conf3) +
                                              (b[[4]]["b_A1_Y"] * sex) +
                                              (b[[4]]["b_A2_Y"] * env) +
                                              (b[[4]]["b_A1A2_Y"] * sex * env)),
                                 sd = b[[5]]["se_Y"])
  }

  data.sim <- data.frame(conf1, conf2, conf3, sex, env, hlth.outcome)
  return(data.sim)
}
