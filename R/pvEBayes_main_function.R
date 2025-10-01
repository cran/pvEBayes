#' Estimate expected null baseline count based on reference row and column
#'
#' @description
#' This function estimates the expected null baseline count (\eqn{E_{ij}}) for
#' each AE-drug combination under the assumption of independence between rows and
#' columns. The expected count is calculated using a reference row (other AEs)
#' and reference column (other drugs). This null baseline is typically used in
#' the empirical Bayes modeling of \code{pvEBayes} package for signal detection
#' and estimation in spontaneous reporting system (SRS) data.
#'
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns). The
#' reference row "Other AEs" and the reference column "Other drugs" need to be the
#' I-th row and J-th column respectively.
#'
#' @details
#' This null value estimator is proposed by Tan et al. (2025).
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches to
#' Pharmacovigilance for Simultaneous Signal Detection and Signal Strength Estimation
#' in Spontaneous Reporting Systems Data. \emph{arXiv preprint.} 2025; arXiv:2502.09816.
#'
#' @return
#'
#' an \code{nrow(contin_table)} by \code{ncol(contin_table)} matrix.
#'
#' @keywords internal
#'
calculate_tilde_e <- function(contin_table) {
  contin_table <- as.matrix(contin_table)
  stopifnot(
    is.numeric(contin_table),
    all(contin_table >= 0),
    all(contin_table == floor(contin_table))
  )
  I <- nrow(contin_table)
  J <- ncol(contin_table)
  res <- contin_table[, J] %*% t(contin_table[I, ]) /
    (contin_table[I, J])
  res %>%
    magrittr::set_rownames(rownames(contin_table)) %>%
    magrittr::set_colnames(colnames(contin_table))
}

.is_valid_contin_table <- function(contin_table) {
  is_valid <- is.numeric(contin_table) &
    all(contin_table >= 0) &
    (all(contin_table == floor(contin_table)))

  if (is_valid == FALSE) {
    warning("contin_table must be a matrix with each entry being non-negative integer")
  }
  return(is_valid)
}



#' Estimate expected null baseline count based on reference row and column
#'
#' @description
#' This function estimates the expected null baseline count (\eqn{E_{ij}}) for
#' each AE-drug combination under the assumption of independence between rows and
#' columns. The expected count is calculated using a reference row (other AEs)
#' and reference column (other drugs). This null baseline is typically used in
#' empirical Bayes modeling of \code{pvEBayes} package for signal detection
#' and estimation in spontaneous reporting system (SRS) data.
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns). The
#' reference row "Other AEs" and the reference column "Other drugs" need to be the
#' I-th row and J-th column respectively.
#'
#' @details
#' This null value estimator is proposed by Tan et al. (2025).
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches to
#' Pharmacovigilance for Simultaneous Signal Detection and Signal Strength Estimation
#' in Spontaneous Reporting Systems Data. \emph{arXiv preprint.} 2025; arXiv:2502.09816.
#'
#' @return
#'
#' an \code{nrow(contin_table)} by \code{ncol(contin_table)} matrix.
#'
#' @export
#'
#' @examples
#'
#'
#' estimate_null_expected_count(statin2025_44)
#'
estimate_null_expected_count <- function(contin_table) {
  calculate_tilde_e(contin_table)
}



.set_default_names <- function(contin_table) {
  I <- nrow(contin_table)
  J <- ncol(contin_table)
  AE_names <- vapply(1:I, function(e) {
    glue::glue("AE", e)
  }, FUN.VALUE = character(1))
  AE_names[I] <- "Other_AEs"

  drug_names <- vapply(1:J, function(e) {
    glue::glue("drug", e)
  }, FUN.VALUE = character(1))
  drug_names[J] <- "Other_drugs"

  contin_table %>%
    magrittr::set_rownames(AE_names) %>%
    magrittr::set_colnames(drug_names)
}

.seq_sobol <- function(from, to, length.out) {
  sobol_seq <- SobolSequence::sobolSequence.points(
    dimR = 2, dimF2 = 14,
    count = length.out
  )[, 1]
  range <- to - from
  res <- sobol_seq * range + from
  return(res)
}


.grid_based_on_hist_log_scale_sobol <- function(N, E, max_draws = FALSE) {
  tmp <- graphics::hist(as.vector(log(N / E)), plot = FALSE)
  num_seq <- length(tmp$breaks) - 1
  num_zero <- sum(N == 0)
  if (max_draws == FALSE) {
    tmp1 <- unique(c(N / E))
    if (length(tmp1) > 1000) {
      return(tmp1)
    }
    c <- 10
  } else {
    c <- max_draws / length(c(N))
  }
  seqs <- lapply(1:num_seq, function(e) {
    .seq_sobol(
      from = tmp$breaks[e], tmp$breaks[e + 1],
      length.out = round(tmp$counts[e] * c)
    )
  }) %>%
    unlist() %>%
    exp()


  if (num_zero != 0) {
    grid_zero <- .seq_sobol(
      from = log(10^(-10)),
      to = tmp$breaks[1],
      length.out = round(num_zero * c)
    )

    seqs <- c(exp(grid_zero), seqs) %>% sort()
  }

  return(seqs)
}

.sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

.get_initial_h <- function(grid) {
  dist21 <- abs(grid - 1)
  weights <- 2 - 2 * .sigmoid(dist21 / 2)
  h <- 1e-7 + (1e-4 - 1e-7) * (weights)
  return(h)
}



.KM_fit <- function(N, E) {
  if (!requireNamespace("Rmosek", quietly = TRUE)) {
    stop("The 'Rmosek' package is required for KM model fitting. Please install it separately.")
  }
  grid <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = 3010)
  fit <- REBayes::Pmix(as.vector(N),
    v = grid, exposure = as.vector(E),
    rtol = 1e-6
  )
  g <- fit$y


  out <- list(g = g, grid = grid, loglik = fit$logLik)
  return(out)
}

.E_fit <- function(N, E, c0 = 1, pDegree = 5, aStart = 1, ...) {
  tau <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = 3010)
  x <- as.vector(N)
  E <- as.vector(E)
  P <- vapply(tau, function(ee) {
    stats::dpois(x, ee * E)
  }, FUN.VALUE = numeric(length(x)))
  m <- length(tau)
  Q <- cbind(1, scale(splines::ns(tau, pDegree),
    center = TRUE,
    scale = FALSE
  ))
  Q <- apply(Q, 2, function(w) w / sqrt(sum(w * w)))




  loglik_1 <- function(a) {
    g0 <- Q %*% a
    g0_max <- max(g0)
    g <- exp(g0 - g0_max)
    # g <- exp(Q %*% a)
    g <- as.vector(g / sum(g))

    f <- as.vector(P %*% g)
    value <- -sum(log(f)) + c0 * sum(a^2)^0.5
    value
  }
  grad_fun <- function(a) {
    g0 <- Q %*% a
    g0_max <- max(g0)
    g <- exp(g0 - g0_max)
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    Pt <- P / f
    W <- g * (t(Pt) - 1)
    qw <- crossprod(Q, W) # t(Q) %*% W
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a / aa
    qw_rowsums <- rowSums(qw)
    grad <- -qw_rowsums + sDot
    grad
  }
  hess_fun <- function(a) {
    g0 <- Q %*% a
    g0_max <- max(g0)
    g <- exp(g0 - g0_max)
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    Pt <- P / f
    W <- g * (t(Pt) - 1)
    qw <- crossprod(Q, W) # t(Q) %*% W
    aa <- sqrt(sum(a^2))
    qw_rowsums <- rowSums(qw)
    qg <- crossprod(Q, g)
    t1 <- tcrossprod(qw)
    t2 <- tcrossprod(qw_rowsums, qg)
    t3 <- t(t2)
    W_rowsums <- rowSums(W)
    t4 <- crossprod(Q * W_rowsums, Q)
    sDotDot <- c0 / aa * (diag(length(a)) - tcrossprod(a) / aa^2)
    hess <- (t1 + t2 + t3 - t4) + sDotDot
    hess
  }
  loglik_not_pena <- function(a) {
    g <- exp(Q %*% a)
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    -sum(log(f))
  }

  aStart <- rep(aStart, ncol(Q))


  result <- stats::nlminb(
    start = aStart, objective = loglik_1, gradient = grad_fun,
    hessian = hess_fun
  )



  mle <- result$par

  g <- exp(Q %*% mle)
  g <- as.vector(g / sum(g))


  hess_fun_non_penal <- function(a) {
    g0 <- Q %*% a
    g0_max <- max(g0)
    g <- exp(g0 - g0_max)
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    Pt <- P / f
    W <- g * (t(Pt) - 1)
    qw <- crossprod(Q, W) # t(Q) %*% W
    aa <- sqrt(sum(a^2))
    qw_rowsums <- rowSums(qw)
    qg <- crossprod(Q, g)
    t1 <- tcrossprod(qw)
    t2 <- tcrossprod(qw_rowsums, qg)
    t3 <- t(t2)
    W_rowsums <- rowSums(W)
    t4 <- crossprod(Q * W_rowsums, Q)
    # sDotDot <- c0/aa * (diag(length(a)) - tcrossprod(a)/aa^2)
    hess <- (t1 + t2 + t3 - t4) #+ sDotDot
    hess
  }

  Hp <- hess_fun(mle)
  H <- hess_fun_non_penal(mle)

  k <- tryCatch(
    {
      (chol2inv(chol(Hp)) %*% H) %>%
        diag() %>%
        sum()
    },
    error = function(e) Inf
  )

  list(
    a = mle, g = g, grid = tau, loglik = -loglik_not_pena(mle),
    df = k
  )
}






#' @useDynLib pvEBayes, .registration = TRUE
#' @importFrom Rcpp evalCpp
.NBmix_EM <- function(N, E, dirichlet = TRUE,
                      alpha = NULL, K = NULL,
                      maxi = NULL,
                      h = NULL,
                      eps = 1e-4) {
  if (dirichlet == FALSE) {
    if (is.null(K)) {
      stop("Error: The number of mixture components K is not given.")
    }
    if (!is.numeric(K) || K <= 1 || floor(K) != K) {
      stop("Error: K must be an integer greater than 1.")
    }
    if (!is.null(alpha)) {
      message("parameter alpha is not needed in GPS/K-gamma model, now running without alpha")
    }
    alpha <- 1
    if (K == 2) {
      grid <- c(1, 2)
    }
    if (K == 3) {
      grid <- c(0, 1, 2)
    }
    if (K > 3) {
      grid <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = K)
    }
  } else {
    if (is.null(alpha)) {
      stop("Error: Dirichlet parameter alpha (0,1) is not given.")
    }
    if (!is.null(K)) {
      message("parameter K is not needed in general-gamma model, now running without K")
    }

    grid <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = 200)
  }
  if (is.null(maxi)) {
    maxi <- 1000
  }
  if (is.null(h)) {
    h <- .get_initial_h(grid)
  }
  result <- NBmix_s1_EM_g_Rcpp(
    N,
    E,
    h,
    grid,
    alpha,
    maxi,
    eps,
    dirichlet
  )
  return(result)
}


.calc_qn_mix_gamma <- function(pvEBayes_obj, N, E, log = FALSE) {
  num_comp <- length(pvEBayes_obj$r)
  if (num_comp == 1) {
    res <- matrix(1, length(N), 1)
    return(res)
  }
  alphas <- pvEBayes_obj$r
  h <- pvEBayes_obj$h
  P <- pvEBayes_obj$omega
  N <- as.vector(N)
  E <- as.vector(E)
  if (log == FALSE) {
    Tau <- vapply(seq_len(length(N)), function(e) {
      tmp <- stats::dnbinom(N[e],
        size = alphas,
        prob = 1 / (1 + E[e] * h), log = TRUE
      ) + log(P)
      tmp <- exp(tmp - max(tmp))
    }, FUN.VALUE = numeric(length(alphas))) %>% t()
    Qn <- Tau / rowSums(Tau)
  } else {
    Tau <- vapply(seq_len(length(N)), function(e) {
      tmp <- stats::dnbinom(N[e],
        size = alphas,
        prob = 1 / (1 + E[e] * h), log = TRUE
      ) + log(P)
      tmp <- tmp - max(tmp)
    }, FUN.VALUE = numeric(length(alphas))) %>% t()
    Qn <- Tau - log(rowSums(exp(Tau)))
  }

  return(Qn)
}

.generate_posterior_gamma_mix <- function(N,
                                          E,
                                          pvEBayes_obj,
                                          nsim = 1000,
                                          ...) {
  prior_shape <- pvEBayes_obj$r
  prior_rate <- 1 / pvEBayes_obj$h
  m <- length(prior_shape)
  I <- nrow(N)
  J <- ncol(N)

  Qn <- .calc_qn_mix_gamma(pvEBayes_obj, N, E)
  Qn <- array(as.vector(Qn), c(I, J, m))
  lambda_sim <- post_draw_gmix_cpp(
    prior_shape,
    prior_rate,
    c(Qn),
    N, E, nsim
  ) %>%
    aperm(c(3, 1, 2))


  dimnames(lambda_sim) <- list(
    NULL, # No names for dim 1
    rownames(N), # Names for dim 2
    colnames(N) # Names for dim 3
  )

  return(lambda_sim)
}


.generate_posterior_grid_based <- function(N, E,
                                           pvEBayes_obj,
                                           nsim = 1000,
                                           ...) {
  grid <- pvEBayes_obj$grid
  esti_prior <- pvEBayes_obj$g
  lambda_sim <- post_draw_discrete_cpp(
    grid,
    esti_prior,
    N, E, nsim
  ) %>%
    aperm(c(3, 1, 2))

  dimnames(lambda_sim) <- list(
    NULL, # No names for dim 1
    rownames(N), # Names for dim 2
    colnames(N) # Names for dim 3
  )

  return(lambda_sim)
}






#' Fit a general-gamma, GPS, K-gamma, KM or efron model for a contingency table.
#'
#' @description
#' This function fits a non-parametric empirical Bayes model to an AE-drug
#' contingency table using one of several empirical Bayes approaches with
#' specified hyperparameter, if is required. Supported models include the
#' "general-gamma", "GPS", "K-gamma", "KM", and "efron".
#'
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns).
#' @param model the model to fit. Available models are "general-gamma",
#' "K-gamma", "GPS", "KM" and "efron". Default to "general-gamma".
#' @param alpha numeric between 0 and 1. The hyperparameter of "general-gamma" model.
#' It is needed if "general-gamma" model is used.
#' @param K integer greater than or equal to 2. It is needed if "K-gamma" model is used.
#' @param p integer greater than or equal to 2. It is needed if "efron" mode is used.
#' @param c0 numeric and greater than 0. It is needed if "efron" mode is used.
#' @param maxi upper limit of iteration for the ECM algorithm.
#' @param n_posterior_draws number of posterior draws for each AE-drug
#' combination.
#' @param eps a tolerance parameter used in the stopping rule of the ECM algorithm.
#'  If the difference in marginal likelihood between two consecutive iterations is
#'  less than eps, the ECM algorithm stops. Default to be 1e-4.
#'
#' @details
#'
#' This function implements the ECM algorithm proposed by Tan et al. (2025),
#' providing a stable and efficient implementation of Gamma-Poisson Shrinker(GPS),
#' K-gamma and "general-gamma" methods for signal estimation and signal detection
#' in Spontaneous Reporting System (SRS) data table.
#'
#' Method "GPS" is proposed by DuMouchel (1999) and it is a parametric empirical
#' Bayes model with a two gamma mixture prior distribution.
#'
#' Methods "K-gamma" and "general-gamma" are non-parametric empirical Bayes models,
#' introduced by Tan et al. (2025). The number of mixture components "K" needs
#' to be prespecified when fitting a "K-gamma" model. For "general-gamma", the
#' mixture weights are regularized by a Dirichlet hyper prior with hyperparameter
#' \eqn{0 \leq \alpha < 1} that controls the shrinkage strength. As "alpha" approaches
#' 0, less non-empty mixture components exist in the fitted model. When \eqn{\alpha = 0},
#' the Dirichlet distribution is an improper prior still offering a reasonable
#' posterior inference that represents the strongest shrinkage of the "general-gamma"
#' model.
#'
#' The implementation of the "KM" model relies on the \pkg{REBayes} package.
#' The model fitting requires the MOSEK optimization solver. Please ensure that
#' \pkg{Rmosek} is correctly installed and configured.
#'
#'
#' The implementation of the "efron" model in this package is adapted from the
#' \pkg{deconvolveR} package, developed by Bradley Efron and
#' Balasubramanian Narasimhan. The original implementation in \pkg{deconvolveR}
#' does not support an exposure or offset parameter in the Poisson model,
#' which corresponds to the expected null value (\eqn{E_{ij}}) for each AE-drug combination.
#' To address this, we modified the relevant code to allow for the inclusion
#' of \eqn{E_{ij}} in the Poisson likelihood. In addition, we implemented a method for
#' estimating the degrees of freedom, enabling AIC- or BIC-based hyperparameter
#' selection for the "efron" model (Tan et al. 2025).
#' See \code{\link{pvEBayes_tune}} for details.
#'
#'
#'
#' @references
#'
#' DuMouchel W. Bayesian data mining in large frequency tables, with an
#' application to the FDA spontaneous reporting system. \emph{The American Statistician.}
#' 1999; 1;53(3):177-90. \cr
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches to
#' Pharmacovigilance for Simultaneous Signal Detection and Signal Strength Estimation
#' in Spontaneous Reporting Systems Data. \emph{arXiv preprint.} 2025; arXiv:2502.09816.
#'
#' Narasimhan B, Efron B. deconvolveR: A G-modeling program for deconvolution
#' and empirical Bayes estimation. \emph{Journal of Statistical Software}.
#' 2020; 2;94:1-20.
#'
#' Koenker R, Gu J. REBayes: an R package for empirical Bayes mixture methods.
#' \emph{Journal of Statistical Software}. 2017; 4;82:1-26.
#'
#'
#'
#'
#'
#'
#' @return
#'
#' The function returns an S3 object of class `pvEBayes` containing the
#' estimated model parameters as well as posterior draws for each AE-drug
#' combination if the number of posterior draws is specified.
#'
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' ref_table <- statin2025_44
#'
#' # set up signal matrix with one signal at entry (1,1)
#' sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
#' sig_mat[1, 1] <- 2
#'
#' # set up structural zero matrix
#' Z <- matrix(0, nrow(ref_table), ncol(ref_table))
#' Z[5, 1] <- 1
#'
#' simu_table <- generate_contin_table(
#'   ref_table = ref_table,
#'   signal_mat = sig_mat,
#'   n_table = 1,
#'   Variation = TRUE,
#'   zi_indic_mat = Z
#' )[[1]][[1]]
#'
#'
#' # fit general-gamma model with a specified alpha
#' fit <- pvEBayes(
#'   contin_table = simu_table, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = 1000
#' )
#'
pvEBayes <- function(contin_table, model = "general-gamma",
                     alpha = NULL, K = NULL,
                     p = NULL, c0 = NULL,
                     maxi = NULL,
                     eps = 1e-4,
                     n_posterior_draws = 1000) {
  h <- NULL
  if (.is_valid_contin_table(contin_table) == FALSE) {
    stop()
  }
  if (is.null(colnames(contin_table)) |
    is.null(rownames(contin_table))) {
    contin_table <- .set_default_names(contin_table)
  }
  E <- calculate_tilde_e(contin_table)
  stopifnot(
    length(model) == 1,
    (model %in% c("general-gamma", "K-gamma", "GPS", "KM", "efron"))
  )
  if (model == "general-gamma") {
    dirichlet <- TRUE
  } else {
    dirichlet <- FALSE
    if (model == "GPS") {
      K <- 2
    }
  }
  start_time <- Sys.time()
  if (model == "KM") {
    res <- .KM_fit(
      contin_table, E
    )
    generate_posterior_fun <- .generate_posterior_grid_based
  } else if (model == "efron") {
    if (is.null(c0) | is.null(p)) {
      stop("Error: hyperparameters c0 and p are not given")
    }
    res <- .E_fit(
      contin_table, E,
      c0, p
    )
    res$c0 <- c0
    res$p <- p
    generate_posterior_fun <- .generate_posterior_grid_based
  } else {
    res <- .NBmix_EM(
      contin_table, E, dirichlet,
      alpha, K,
      maxi,
      h,
      eps
    )
    res$alpha <- alpha
    generate_posterior_fun <- .generate_posterior_gamma_mix
  }

  end_time <- Sys.time() # start the clock

  res$fit_time <- difftime(end_time, start_time)
  res$draws_time <- NULL
  res$model <- model

  res$n_posterior_draws <- n_posterior_draws
  if (!is.null(n_posterior_draws)) {
    stopifnot(
      is.numeric(n_posterior_draws),
      n_posterior_draws > 0,
      n_posterior_draws == round(n_posterior_draws)
    )

    res$posterior_draws <- generate_posterior_fun(contin_table,
      E,
      res,
      nsim = n_posterior_draws
    )
    end_time2 <- Sys.time()
    res$draws_time <- difftime(end_time2, end_time)
  }
  res$contin_table <- contin_table
  res$E <- E

  class(res) <- "pvEBayes"
  res
}


#' Select hyperparameter and obtain the optimal general-gamma or efron model
#'  based on AIC and BIC
#'
#' @description
#' This function performs hyperparameter tuning for the general-gamma or Efron
#' model using AIC or BIC. For a given AE-drug contingency table, the
#' function fits a series of models across a grid of candidate hyperparameter values
#' and computes their AIC and BIC. The models with the lowest AIC or BIC values
#' are selected as the optimal fits.
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns).
#' @param model the model to be tuned. Available models are "general-gamma" and
#' "efron". Default to "general-gamma".
#' @param alpha_vec vector of hyperparameter alpha values to be selected. Alpha
#' is a hyperparameter in general-gamma model which is numeric and between 0 and 1.
#' If is NULL, a default set of alpha values (0, 0.1, 0.3, 0.5, 0.7, 0.9) will
#' be used.
#' @param p_vec vector of hyperparameter p values to be selected. p is
#' a hyperparameter in "efron" model which should be a positive integer. If is NULL,
#' a default set of p values (80, 100, 120, 150, 200) will be used.
#' @param c0_vec vector of hyperparameter c0 values to be selected. c0 is
#' a hyperparameter in "efron" model which should be a positive number. If is NULL,
#' a default set of c0 values (0.001, 0.01, 0.1, 0.2, 0.5) will be used.
#' @param use_AIC logical, indicating whether AIC or BIC is used. Default to be
#' TRUE.
#' @param n_posterior_draws number of posterior draws for each AE-drug
#' combination.
#' @param return_all_fit logical, indicating whether to return all the fitted
#' model under the selection. Default to be FALSE.
#' @param return_all_AIC logical, indicating whether to return AIC values for
#' each fitted model under the selection. Default to be TRUE.
#' @param return_all_BIC logical, indicating whether to return BIC values for
#' each fitted model under the selection. Default to be TRUE.
#'
#'
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control.}
#' 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model. \emph{The Annals of Statistics.}
#' 1978; 1:461-4.
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches to
#' Pharmacovigilance for Simultaneous Signal Detection and Signal Strength Estimation
#' in Spontaneous Reporting Systems Data. \emph{arXiv preprint.} 2025; arXiv:2502.09816.
#'
#' @return
#'
#' The function returns an S3 object of class `pvEBayes` containing the selected
#' estimated model parameters as well as posterior draws for each AE-drug
#' combination if the number of posterior draws is specified.
#'
#' @export
#'
#' @examples
#'
#' fit <- pvEBayes_tune(statin2025_44, model = "general-gamma")
#'
pvEBayes_tune <- function(contin_table, model = "general-gamma",
                          alpha_vec = NULL,
                          p_vec = NULL, c0_vec = NULL,
                          use_AIC = TRUE,
                          n_posterior_draws = 1000,
                          return_all_fit = FALSE,
                          return_all_AIC = TRUE,
                          return_all_BIC = TRUE) {
  if (.is_valid_contin_table(contin_table) == FALSE) {
    stop()
  }
  if (is.null(colnames(contin_table)) |
    is.null(rownames(contin_table))) {
    contin_table <- .set_default_names(contin_table)
  }
  E <- calculate_tilde_e(contin_table)
  stopifnot(
    length(model) == 1,
    (model %in% c("general-gamma", "K-gamma", "GPS", "KM", "efron"))
  )

  if (model %in% c("K-gamma", "GPS", "KM")) {
    stop("Please use pvEBayes() for GPS, K-gamma or KM model fitting.")
  }

  if (model == "efron") {
    objects <- tuning_efron(
      contin_table = contin_table,
      p_vec = p_vec,
      c0_vec = c0_vec,
      return_all_fit = return_all_fit,
      return_all_AIC = return_all_AIC,
      return_all_BIC = return_all_BIC
    )
    generate_posterior_fun <- .generate_posterior_grid_based
  } else {
    objects <- tuning_general_gamma(
      contin_table = contin_table,
      alpha_vec = alpha_vec,
      return_all_fit = return_all_fit,
      return_all_AIC = return_all_AIC,
      return_all_BIC = return_all_BIC
    )
    generate_posterior_fun <- .generate_posterior_gamma_mix
  }


  if (use_AIC == TRUE) {
    res <- objects$best_model_AIC
  } else {
    res <- objects$best_model_BIC
  }

  end_time <- Sys.time() # END the clock

  res$draws_time <- NULL
  res$model <- model

  res$n_posterior_draws <- n_posterior_draws
  if (!is.null(n_posterior_draws)) {
    stopifnot(
      is.numeric(n_posterior_draws),
      n_posterior_draws > 0,
      n_posterior_draws == round(n_posterior_draws)
    )

    res$posterior_draws <- generate_posterior_fun(contin_table,
      E,
      res,
      nsim = n_posterior_draws
    )
    end_time2 <- Sys.time()
    res$draws_time <- difftime(end_time2, end_time)
  }
  res$contin_table <- contin_table
  res$E <- E

  res$tuning <- objects

  class(res) <- c("pvEBayes", "pvEBayes_tuned")
  res
}




#' Select hyperparameter alpha and obtain the optimal general-gamma model based
#'  on AIC and BIC
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns).
#' @param alpha_vec vector of hyperparameter alpha values to be selected. Alpha
#' is hyperparameter in general-gamma model which is numeric and between 0 and 1.
#' If is NULL, a default set of alpha values (0, 0.1, 0.3, 0.5, 0.7, 0.9) will
#' be used.
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control.}
#' 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model. \emph{The Annals of Statistics.}
#' 1978; 1:461-4.
#'
#'
#' @return
#' a list of fitted models with hyperparameter alpha selected by AIC or BIC.
#'
#' @keywords internal
#'
tuning_general_gamma <- function(contin_table,
                                 alpha_vec = NULL,
                                 return_all_fit = FALSE,
                                 return_all_AIC = TRUE,
                                 return_all_BIC = TRUE) {
  if (is.null(alpha_vec)) {
    alpha_vec <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9)
  }
  stopifnot(
    length(alpha_vec) > 0,
    is.numeric(alpha_vec),
    all(alpha_vec >= 0 & alpha_vec <= 1),
    .is_valid_contin_table(contin_table)
  )
  fits <- alpha_vec %>%
    lapply(function(e) {
      pvEBayes(
        contin_table = contin_table, model = "general-gamma",
        alpha = e, n_posterior_draws = NULL,
        eps = 1e-4
      )
    })
  AICs <- fits %>% vapply(AIC.pvEBayes, FUN.VALUE = numeric(1))
  BICs <- fits %>% vapply(BIC.pvEBayes, FUN.VALUE = numeric(1))
  num_mixture <- fits %>% vapply(function(e) {
    length(e$r)
  }, FUN.VALUE = numeric(1))

  res <- data.frame(
    alpha = alpha_vec,
    AIC = AICs,
    BIC = BICs,
    num_mixture = num_mixture
  )

  AIC_winner <- res$alpha[which.min(res$AIC)]
  BIC_winner <- res$alpha[which.min(res$BIC)]
  message(
    glue::glue(
      "The alpha value selected under AIC is ", AIC_winner, ",\n",
      "The alpha value selected under BIC is ", BIC_winner, "."
    ),
    "\n"
  )
  message(paste(utils::capture.output(res), collapse = "\n"))

  out <- list(
    best_model_AIC = fits[[which.min(res$AIC)]],
    best_model_BIC = fits[[which.min(res$BIC)]]
  )

  if (return_all_fit) {
    out$all_fit <- fits
  }

  if (return_all_AIC) {
    out$all_AIC <- AICs
  }

  if (return_all_BIC) {
    out$all_BIC <- BICs
  }

  invisible(out)
}


.get_grid_list <- function(c0, p) {
  res <- list()
  for (i in seq_len(length(c0))) {
    for (j in seq_len(length(p))) {
      res[[10 * i + j]] <- list(c0 = c0[i], p = p[j])
    }
  }
  res <- Filter(Negate(is.null), res)
  return(res)
}


#' Select hyperparameter (p, c0) and obtain the optimal efron model based
#'  on AIC and BIC
#'
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of adverse
#' events for I AEs (along the rows) and J drugs (along the columns).
#' @param p_vec vector of hyperparameter p values to be selected. p is
#' a hyperparameter in "efron" model which should be a positive integer. If is NULL,
#' a default set of p values (80, 100, 120, 150, 200) will be used.
#' @param c0_vec vector of hyperparameter c0 values to be selected. c0 is
#' a hyperparameter in "efron" model which should be a positive number. If is NULL,
#' a default set of c0 values (0.001, 0.01, 0.1, 0.2, 0.5) will be used.
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control.}
#' 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model. \emph{The Annals of Statistics.}
#' 1978; 1:461-4.
#'
#' @return
#' a list of fitted models with hyperparameter alpha selected by AIC or BIC.
#'
#' @keywords internal
#'
tuning_efron <- function(contin_table,
                         p_vec = NULL,
                         c0_vec = NULL,
                         return_all_fit = FALSE,
                         return_all_AIC = TRUE,
                         return_all_BIC = TRUE) {
  if (is.null(p_vec)) {
    p_vec <- c(40, 60, 80, 100, 120)
  }
  if (is.null(c0_vec)) {
    c0_vec <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
  }
  stopifnot(
    length(p_vec) > 0,
    length(c0_vec) > 0,
    all(p_vec >= 2),
    all(c0_vec >= 0),
    .is_valid_contin_table(contin_table)
  )
  parameter_list <- .get_grid_list(c0_vec, p_vec)
  fits <- parameter_list %>%
    lapply(function(e) {
      pvEBayes(
        contin_table = contin_table, model = "efron",
        p = e$p, c0 = e$c0,
        n_posterior_draws = NULL
      )
    })
  AICs <- fits %>% vapply(AIC.pvEBayes, FUN.VALUE = numeric(1))
  BICs <- fits %>% vapply(BIC.pvEBayes, FUN.VALUE = numeric(1))


  res <- data.frame(
    p = parameter_list %>% vapply(function(e) e$p, FUN.VALUE = numeric(1)),
    c0 = parameter_list %>% vapply(function(e) e$c0, FUN.VALUE = numeric(1)),
    AIC = AICs,
    BIC = BICs
  )

  AIC_winner <- (res[, 1:2][which.min(res$AIC), ]) %>%
    unlist() %>%
    unname()
  BIC_winner <- (res[, 1:2][which.min(res$BIC), ]) %>%
    unlist() %>%
    unname()
  message(
    glue::glue(
      "The hyperparameters selected under AIC is (p = ", AIC_winner[1],
      ", c0 = ", AIC_winner[2], ")", ",\n",
      "The hyperparameters selected under BIC is (p = ", BIC_winner[1],
      ", c0 = ", BIC_winner[2], ").", ",\n",
    ),
    "\n"
  )
  message(paste(utils::capture.output(res), collapse = "\n"))

  out <- list(
    best_model_AIC = fits[[which.min(res$AIC)]],
    best_model_BIC = fits[[which.min(res$BIC)]]
  )

  if (return_all_fit) {
    out$all_fit <- fits
  }

  if (return_all_AIC) {
    out$all_AIC <- AICs
  }

  if (return_all_BIC) {
    out$all_BIC <- BICs
  }
  invisible(out)
}
