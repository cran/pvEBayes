#' NA_standards
#'
#' @srrstatsVerbose TRUE
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @srrstatsNA {G2.4d} Factor input is not used.
#' @srrstatsNA {G2.4e} Factor input is not used.
#' @srrstatsNA {G2.5} Factor input is not used.
#' @srrstatsNA {G2.11} No data.frame-like input is used.
#' @srrstatsNA {G2.12} No data.frame-like input is used.
#'
#'
#' @srrstatsNA {G3.1} There is no covariance computation.
#' @srrstatsNA {G3.1a} There is no covariance computation.
#' @srrstatsNA {G4.0} There is not function that enables outputs to be written
#' to local files.
#' @srrstatsNA {G5.4c} This is not applicable.
#' @srrstatsNA {G5.10, G5.11, G5.11a, G5.12}
#' These tests are not applicable to \code{pvEBayes}. Therefore, are not
#' included.
#'
#' @srrstatsNA {BS1.3, BS1.3a, BS1.3b}
#' The methods implemented in \pkg{pvEBayes} are empirical Bayes methods that
#' are fitted by estimating prior distributions through optimization with
#' respect to joint marginal likelihood rather than
#' through sampling procedures. As a result, the package does not involve
#' stochastic sampling processes and thinning parameters, burn-in lengths,
#' random seeds, or convergence diagnostics for full Bayesian computations
#' are not used. Parameters that control the empirical Bayes optimization are
#' appropriately documented.
#'
#' @srrstatsNA {BS1.5}
#' There is no multiple convergence criteria implemented for each optimizer.
#'
#' @srrstatsNA {BS2.2, BS2.3, BS2.4, BS2.5}
#' \pkg{pvEBayes} is based on parametric/nonparametric empirical Bayes methods
#' in which the prior distributions are estimated from data, not specified by
#' the user. Therefore, users do not provide prior distributional parameters in
#' the sense of full Bayesian modeling.
#'
#' @srrstatsNA {BS2.7, BS2.8, BS2.9, BS2.10, BS2.11}
#' \pkg{pvEBayes} is based on parametric/nonparametric empirical Bayes methods
#' that the methods are fitted by estimating prior distributions through
#' optimization rather than through sampling procedures. Tan et al.
#' (Stat. in Med., 2025) shows that the non-parametric empirical Bayes methods
#' included in \pkg{pvEBayes} are not as computational intensive as the
#' non-parametric full Baysian competitor. Parallel computation is not used.
#'
#' In \pkg{pvEBayes}, initialization of prior parameters follows the
#' recommendations of Tan et al. (*Statistics in Medicine*, 2025).
#' Users are not required to provide initial values for prior parameters,
#' as the package automatically performs this step internally.
#'
#'
#' @srrstatsNA {BS3.1, BS3.2}
#' The models implemented in \pkg{pvEBayes} are empirical Bayes models applied
#' to SRS contingency tables and do not involve regression predictors or
#' design matrices. As such, the notion of collinearity is not applicable in
#' this context.
#'
#' @srrstatsNA {BS4.0, BS4.1, BS4.2, BS4.3, BS4.4, BS4.5, BS4.6, BS4.7}
#' The methods implemented in \pkg{pvEBayes} are parametric/nonparametric
#' empirical Bayes approaches rely on optimizing the joint marginal likelihood
#' with optimizers. These models do not rely on Monte Carlo sampling in full
#' Bayesian modeling.
#'
#' @srrstatsNA {BS5.0, BS5.2}
#' The methods implemented in \pkg{pvEBayes} are parametric/nonparametric
#' empirical Bayes approaches rely on optimizing the joint marginal likelihood
#' with optimizers. Returning values with starting value(s) or seed(s) is not
#' meaningful.
#'
#' @srrstatsNA {BS6.2, BS6.5}
#' The empirical Bayes methods implemented in \pkg{pvEBayes} do not rely on
#' stochastic sampling. Therefore, plotting sequences of posterior samples, with
#' burn-in periods are not enabled.
#'
#' @srrstatsNA {BS7.0, BS7.1, BS7.2}
#' These requirements relate to Bayesian software in which a parametric prior
#' is specified by the user. The methods implemented in \code{pvEBayes} are
#' parametric/non-parametric empirical Bayes approaches in which the prior
#' distributions are estimated from the data.
#'
#' @srrstatsNA {BS7.3} This test is not included in \code{pvEBayes}. The effect
#' of data size on the scaling of algorithmic efficiency is investigated in
#' Tan et al. (Stat. in Med., 2025)
#'
#' @noRd
NULL
