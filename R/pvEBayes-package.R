#' @name pvEBayes-package
#' @title A suite of empirical Bayes methods to use in pharmacovigilance.
#'
#' @description
#' \code{pvEBayes} provides a collection of parametric and non-parametric
#' empirical Bayes methods implementation for pharmacovigilance (including
#' signal detection and signal estimation) on spontaneous reporting systems
#' (SRS) data.
#'
#' An SRS dataset catalogs AE reports on *I* AE rows across *J* drug columns.
#' Let \eqn{N_{ij}} denote the number of reported cases for the
#' *i*-th AE and the *j*-th drug, where \eqn{i = 1,..., I} and
#' \eqn{j = 1,..., J}. We assume that for each AE-drug pair,
#' \eqn{N_{ij} \sim \text{Poisson}(\lambda_{ij} E_{ij})}, where \eqn{E_{ij}} is
#' expected baseline value measuring the expected count of the AE-drug pair
#' when there is no association between *i*-th AE and *j*-th drug. The parameter
#' \eqn{\lambda_{ij} \geq 0} represents the relative reporting ratio, the signal
#' strength, for the \eqn{(i, j)}-th pair measuring the ratio of the actual
#' expected count arising due to dependence to the null baseline expected count.
#' Current disproportionality analysis mainly focuses on \emph{signal detection}
#' which seeks to determine whether the observation \eqn{N_{ij}} is substantially
#' greater than the corresponding null baseline \eqn{E_{ij}}. Under the Poisson
#' model, that is to say, its signal strength \eqn{\lambda_{ij}} is
#' significantly greater than 1.
#'
#' In addition to *signal detection*, Tan et al. (*Stat. in Med.*, 2025) broaden
#' the role of disproportionality to *signal estimation*. The use of the flexible
#' non-parametric empirical Bayes models enables more nuanced empirical Bayes
#' posterior inference (parameter estimation and uncertainty quantification) on
#' signal strength parameter \eqn{\{ \lambda_{ij} \}}. This allows researchers to
#' distinguish AE-drug pairs that would appear similar under a binary signal
#' detection framework. For example, the AE-drug pairs with signal strengths of
#' 1.5 and 4.0 could both be significantly greater than 1 and detected as a
#' signal. Such differences in signal strength may have distinct implications in
#' medical and clinical contexts.
#'
#' The methods included in \code{pvEBayes} differ by their assumptions on the
#' prior distribution. Implemented methods include the Gamma-Poisson Shrinker
#' (GPS), Koenker-Mizera (KM) method, Efron’s nonparametric empirical Bayes
#' approach, the K-gamma model, and the general-gamma model.
#'
#' The GPS model uses two gamma mixture prior by assuming the signal/non-signal
#' structure in SRS data. However, in real-world setting, signal
#' strengths \eqn{(\lambda_{ij})} are often heterogeneous and thus follows a
#' multi-modal distribution, making it difficult to assume a parametric prior.
#' Non-parametric empirical Bayes models (KM, Efron, K-gamma and general-gamma)
#' address this challenge by utilizing a flexible prior with general mixture
#' form and estimating the prior distribution in a data-driven way.
#'
#'
#' \code{pvEBayes} offers the first implemention of the bi-level Expectation
#' Conditional Maximization (ECM) algorithm proposed by Tan et al. (2025) for
#' efficient parameter estimation in gamma mixture prior based models: GPS
#' K-gamma and general-gamma.
#'
#' The KM method has an existing implementation in the \code{REBayes} package,
#' but it relies on Mosek, a commercial convex optimization solver, which may
#' limit accessibility due to licensing issue. \code{pvEBayes} provides a
#' alternative fully open-source implementation of the KM method using
#' \code{CVXR}.
#'
#' Efron’s method also has a general nonparametric empirical Bayes
#' implementation in the \code{deconvolveR} package; however, that
#' implementation does not support an exposure or offset parameter in the
#' Poisson model, which corresponds to the expected null value \eqn{E_{ij}}.
#' In \code{pvEBayes}, the implementation of the Efron's method is adapted and
#' modified from \code{deconvolveR} to support \eqn{E_{ij}} in Poisson model.
#'
#'
#' For a detailed introduction to \code{pvEBayes}, see Tan et al.
#' (\emph{arxiv}, 2025) and package Vignette.
#'
#' @author
#' Yihao Tan, Marianthi Markatou, Saptarshi Chakraborty and Raktim Mukhopadhyay.
#'
#' Maintainer: Yihao Tan \email{yihaotan@buffalo.edu}
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
#'
#' Tan Y, Markatou M and Chakraborty S. pvEBayes: An R Package for Empirical
#' Bayes Methods in Pharmacovigilance.
#' \emph{arXiv}:2512.01057 (stat.AP). https://doi.org/10.48550/arXiv.2512.01057
#'
#' Koenker R, Mizera I. Convex Optimization, Shape Constraints, Compound
#' Decisions, and Empirical Bayes Rules. \emph{Journal of the American
#' Statistical Association} 2014; 109(506): 674–685,
#' https://doi.org/10.1080/01621459.2013.869224
#'
#' Efron B. Empirical Bayes Deconvolution Estimates. \emph{Biometrika} 2016;
#' 103(1); 1-20, https://doi.org/10.1093/biomet/asv068
#'
#' DuMouchel W. Bayesian data mining in large frequency tables, with an
#' application to the FDA spontaneous reporting system.
#' \emph{The American Statistician.} 1999; 1;53(3):177-90.
#'
#' @importFrom magrittr %>%
#' @importFrom data.table .SD
#'
#' @srrstats {G1.0} Reference section reports the related literature.
#' @srrstats {G1.1} The bi-level Expectation Conditional Maximization (ECM)
#' algorithm for gamma mixture based prior models (including GPS, K-gamma and
#' general-gamma) It is the first implementation of a novel algorithm. Our
#' implementation for Koenker-Mizera (KM) and Efron methods are improvement
#' to \code{REBayes} and \code{deconvolveR} package respectively to fit into
#' disproportionality analysis on spontaneous reporting systems SRS data.
#' @srrstats {G1.2} see the 'CONTRIBUTING.md' for states of development.
#' @srrstats {G1.3} All statistical terminologies are clarified and
#' unambiguously defined
#' @srrstats {G1.4} All functions are documented with roxygen2.
#' @srrstats {G1.4a} All internal functions are documented with roxygen2.
#' @srrstats {G1.5} Vignette reproduces the results in the associated
#' publications.
#' @srrstats {G5.0} Several FAERS datasets are included and they are
#' illustrated in example and vignette.
#' @srrstats {G5.1} Relevant datasets are provided. See 'data.R'.
#' @srrstats {G5.2, G5.2a, G5.2b} All exported functions have sufficient checks
#' for inputs and appropriate error/warning messages after that.
#' @srrstats {G5.3} Functions that could potentially return objects containing
#' (`NA`) or undefined (`NaN`, `Inf`) values are tested in
#' 'test-pvEBayes_main_function.R'.
#' @srrstats {G5.8, G5.8a, G5.8b, G5.8c, G5.8d}
#' Edge condition tests are provided in 'test-pvEBayes_main_function.R'.
#' @srrstats {BS1.0} The meaning and effect of the hyperparameters in K-gamma,
#' general-gamma and Efron methods are carefully explained in the function
#' documentation, README and vignette.
#' @srrstats {BS1.1} Descriptions of how to enter data are presented in both
#' textual and code form in documentation and vignette.
#' @srrstats {BS1.2, BS1.2a, BS1.2b, BS1.2c}
#' The methods implemented in \pkg{pvEBayes} use nonparametric empirical Bayes
#' priors whose specification is controlled through model choices and
#' hyperparameters (e.g., \code{alpha}, \code{c0}, \code{K}). A clear
#' description of how these prior distributions are specified is
#' provided at three levels: (i) example in the package README, (ii) examples
#' in the package vignette, and (iii) function-level documentation.
#' @srrstats {BS1.4}
#' Fitting models in \pkg{pvEBayes} do not involve stochastic sampling
#' procedures. Consequently, the package does not implement convergence checkers
#' that are commonly used to assess the convergence of a full Bayesian sampler.
#'
#' Estimation in \pkg{pvEBayes} is performed through deterministic optimization
#' methods: convex optimizer for the KM, log marginal likelihood optimizer for
#' Efron model, and an ECM algorithm for the general-gamma model. These
#' procedures have termination criteria determined by optimization tolerances
#' rather than stochastic convergence diagnostics. The convergence criteria
#' rtol_ecm, rtol_KM and rtol_efron are appropriately documented in the
#' function-level documentation.
#' @srrstats {BS2.1, BS2.1a} \pkg{pvEBayes} implements appropriate input checks
#' to ensure that all user-supplied data objects are valid. The behavior of
#' these checks are covered by the package’s test suite.
#' @srrstats {BS2.2} In \pkg{pvEBayes}, all hyperparameters are validated and
#' pre-processed before they are submitted to analytic routines.
#' @srrstats {BS2.6} \pkg{pvEBayes} performs explicit checks on the estimated
#' prior parameters and applies corrective steps in each iteration of the ECM
#' algorithm.
#'
#' For the KM and Efron methods, \pkg{pvEBayes} relies on the underlying solvers
#' (CVXR::solve and stats::nlminb). These optimizers do not permit
#' user-level intervention during iteration; invalid parameter values produced
#' by the solver cause NaN or infinite objective values, which leads to the
#' optimizer to terminate with an error.
#'
#' @srrstats {BS2.12, BS2.13}
#' The only function in \pkg{pvEBayes} that produces console output is
#' \code{pvEBayes_tune()}, which prints AIC and BIC values for
#' hyperparameters under selection. This output can be suppressed using
#' \code{suppressMessages()} if desired. The suppression behavior is explicitly
#' tested in \code{test-pvEBayes_main_function.R}.
#'
#'
#' @srrstats {G2.13, G2.14, G2.14a, G2.14b, G2.14c, G2.15, G2.16}
#' Missing value in the context of SRS data is underreporting
#' or zero observation for a A (SRS contingency entry), rather than
#' an \code{NA}, \code{Inf}, or other non-numeric placeholder.
#' Such observations are expected and fully supported by the methods implemented
#' in \pkg{pvEBayes}.
#'
#' @srrstats {BS2.14}
#' The \pkg{pvEBayes} package does not currently give warnings during normal
#' operation. Invalid inputs result in errors via \code{stop()}.
#'
#' @srrstats {BS2.15}
#' Functions in \pkg{pvEBayes} report failures using standard R error
#' conditions via \code{stop()}. These errors can be programmatically captured
#' and processed by users using \code{tryCatch()} or similar mechanisms. This
#' is tested in \code{test-pvEBayes_main_function.R}.
#'
#' @srrstats {BS3.0} Missing value in the context of SRS data is underreporting
#' or zero observation for an SRS contingency entry, rather than
#' an \code{NA}, \code{Inf}, or other non-numeric placeholder.
#' Such observations are expected and fully supported by the methods implemented
#' in \pkg{pvEBayes}.
#'
#' @srrstats {BS5.1} Return values include appropriate metadata on
#' types (or classes) and dimensions of input data.
#'
#' @srrstats {BS5.3, BS5.4, bs5.5}
#' The empirical Bayes methods implemented in \pkg{pvEBayes} do not rely on
#' stochastic sampling, and therefore do not produce the types of
#' convergence diagnostics typically associated with full Bayesian modeling.
#' Convergence in the ECM algorithm is reached (at least to a sub-optimal).
#' This is ensured by monotonically increased log joint marginal likelihood,
#' as proved by Tan et al. (*Stat. in Med*, 2025).
#'
#' @srrstats {BS6.0} \code{print()} function is implemented for \code{pvEBayes}
#' and \code{pvEBayes_tuned} objects.
#'
#' @srrstats {BS6.1} \code{plot()} function is implemented for \code{pvEBayes}
#' and \code{pvEBayes_tuned} objects.
#'
#' @srrstats {BS6.3} \code{eyeplot_pvEBayes()} function implemented for
#' \code{pvEBayes} and \code{pvEBayes_tuned} objects provide a visualization
#' for the estimated posterior distribution for the parameter of interest
#' \eqn{\lambda}.
#'
#' @srrstats {BS6.4} \code{summary()} function is implemented for
#' \code{pvEBayes} and \code{pvEBayes_tuned} objects.
#'
#'
#'
"_PACKAGE"
