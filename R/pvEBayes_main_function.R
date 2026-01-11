#' Estimate expected null baseline count based on reference row and column
#'
#' @description
#' This function estimates the expected null baseline count (\eqn{E_{ij}}) for
#' each AE-drug combination under the assumption of independence between rows
#' and columns. The expected count is calculated using a reference row
#' (other AEs) and reference column (other drugs). This null baseline is
#' typically used in the empirical Bayes modeling of \pkg{pvEBayes} package
#' for signal detection and estimation in spontaneous reporting system (SRS)
#' data.
#'
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' The reference row "Other AEs" and the reference column "Other drugs" need to
#' be the I-th row and J-th column respectively.
#'
#' @details
#' This null value estimator is proposed by Tan et al. (2025).
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
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

#' Check if a input is a valid SRS contingency table
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#'
#' @returns logical.
#' @keywords internal
#' @noRd
.is_valid_contin_table <- function(contin_table) {
  is_valid <- is.numeric(contin_table) &
    all(contin_table >= 0) &
    (all(contin_table == floor(contin_table)))

  if (is_valid == FALSE) {
    warning(
      paste0(
        "contin_table must be a matrix with",
        " each entry being non-negative integer."
      )
    )
  }
  is_valid
}



#' Estimate expected null baseline count based on reference row and column
#'
#' @description
#' This function estimates the expected null baseline count (\eqn{E_{ij}}) for
#' each AE-drug combination under the assumption of independence between rows
#' and columns. The expected count is calculated using a reference row
#' (other AEs) and reference column (other drugs). This null baseline is
#' typically used in empirical Bayes modeling of \pkg{pvEBayes} package for
#' signal detection and estimation in spontaneous reporting system (SRS) data.
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' The reference row "Other AEs" and the reference column "Other drugs" need to
#' be the I-th row and J-th column respectively.
#'
#' @details
#' This null value estimator is proposed by Tan et al. (2025).
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
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



#' Set default row and column names to a contingency table
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#'
#' @returns the contingency table with default row names (AE1, AE2, ...)and
#' column names (drug1, drug2, ...).
#'
#'
#' @keywords internal
#' @noRd
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



#' Generate one-dimensional uniform quasi-random sequence
#'
#' @param from the lower bound of the sequence
#' @param to the upper bound of the sequence
#' @param length.out desired length of the sequence
#'
#' @returns a numeric vector
#' @keywords internal
#' @noRd
.seq_sobol <- function(from, to, length.out) {
  sobol_seq <- numeric(length.out)
  for (i in 0:(length.out-1)) {
    x <- 0
    k <- i
    base <- 0.5
    while (k > 0) {
      if (k %% 2 == 1) x <- x + base
      k <- k %/% 2 # Integer division
      base <- base / 2
    }
    sobol_seq[i] <- x
  }
  range <- to - from
  res <- sobol_seq * range + from
  res
}


#' Histogram-based grid value generation for discrete non-parametric empirical
#' Bayes methods
#'
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E an IxJ numeric table of estimated null baseline count. E should be
#' the same dimension as N.
#' @param max_draws a upper limit of the generated grid size. If a max_draws
#' is provided, the grid size is min(max_draws, 10*I*J). Defaul to NULL.
#'
#' @returns a numeric vector
#' @keywords internal
#' @noRd
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

  seqs
}

#' Sigmoid function
#'
#' @param x numeric input
#'
#' @returns numeric sigmoid return
#' @keywords internal
#' @noRd
.sigmoid <- function(x) {
  1 / (1 + exp(-x))
}


#' Obtain initialization of h for general-gamma method
#'
#' @param grid a numeric vector
#'
#' @returns a numeric vector
#' @keywords internal
#' @noRd
.get_initial_h <- function(grid) {
  dist21 <- abs(grid - 1)
  weights <- 2 - 2 * .sigmoid(dist21 / 2)
  h <- 1e-7 + (1e-4 - 1e-7) * (weights)
  h
}

#' Parameter estimation for Koenker-Mizera method
#'
#' @param x a vector of sample observation. In the context of SRS mining,
#' x is vectorized SRS contingency table.
#' @param v a vector of grid values
#' @param exposure a vector of observation specific exposures. In the context
#' of SRS mining, it is a vector of expected null baseline values.
#' @param ... other parameters passed to the CVXR optimizer.
#'
#' @returns a list of CVXR optimizer outputs
#' @keywords internal
#' @noRd
.km_eb_fit <- function(x, v, exposure = NULL, ...) {
  n <- length(x)

  m <- length(v)
  d <- rep(1, length(v))
  w <- rep(1, n) / n
  A <- matrix(0, n, m)

  for (i in seq_len(n)) {
    A[i, ] <- vapply(
      X = v * exposure[i],
      FUN = function(lambda) {
        stats::dpois(x[i], lambda)
      },
      FUN.VALUE = numeric(1L)
    )
  }

  f <- .KWDual_CVXR(A, d, w, ...)
  logLik <- n * sum(w * log(f$g))
  z <- list(
    x = v, y = f$f, g = f$g, logLik = logLik,
    status = f$status
  )
  z
}



#' Convex optimization for Koenker-Mizera method with CVXR
#'
#' @param A linear constraint matrix. In the context of SRS mining, entries in
#' A represent the conditional Poisson likelihood of the observation given
#' grid value.
#' @param d constraint vector. In the context of SRS mining, d is the
#' reciprocal of grid values.
#' @param w weights for observations.
#' @param rtol_KM relative tolerance for optimization algorithm.
#' @param verb
#'
#' @returns a list of CVXR optimizer outputs
#' @keywords internal
#' @noRd
.KWDual_CVXR <- function(A, d, w, rtol_KM = 1e-6, verb = FALSE) {
  m <- ncol(A)

  A_mat <- A
  fvar <- CVXR::Variable(m)
  gexpr <- A_mat %*% (fvar * d)

  objective <- CVXR::Maximize(CVXR::sum_entries(w * log(gexpr)))
  constraints <- list(
    fvar >= 0,
    CVXR::sum_entries(d * fvar) == 1
  )

  prob <- CVXR::Problem(objective, constraints)
  res <- CVXR::solve(prob,
    solver = "ECOS",
    reltol = rtol_KM,
    verbose = verb
  )

  fhat <- as.vector(res$getValue(fvar))
  fhat[fhat < 0] <- 0
  ghat <- as.vector(A_mat %*% (fhat * d))

  list(f = fhat, g = ghat, status = res$status)
}


#' Fit a Koenker-Mizera (KM) model for a contingency table.
#'
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param rtol_KM The relative tolerance on the duality gap.
#'
#' @details
#'
#' Parameter estimation for the "KM" model is formulated as a convex
#' optimization problem. The objective function and constraints used in
#' \pkg{pvEBayes} follow the same construction as in \pkg{REBayes}.
#' Parameter estimation is performed using the open-source convex optimization
#' package \pkg{CVXR}. The grid value generation follows the recommendation of
#' Tan et al. (2025).
#'
#' @references
#'
#' Koenker R, Gu J. REBayes: an R package for empirical Bayes mixture methods.
#' \emph{Journal of Statistical Software}. 2017; 4;82:1-26.
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
#'
#'
#' Fu, A, Narasimhan, B, Boyd, S. CVXR: An R Package for Disciplined Convex
#' Optimization. \emph{Journal of Statistical Software}. 2020; 94;14:1-34.
#'
#'
#' @returns a list of CVXR optimizer outputs
#' @keywords internal
.KM_fit <- function(N, E, rtol_KM = 1e-6) {
  n_draws <- prod(dim(N)) * 3
  if (n_draws >= 1000) {
    n_draws <- 1000
  }
  grid <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = n_draws)
  fit <- .km_eb_fit(as.vector(N),
    v = grid, exposure = as.vector(E),
    rtol_KM = rtol_KM
  )
  g <- fit$y


  out <- list(g = g, grid = grid, loglik = fit$logLik)
  out
}

#' Fit an Efron model for a contingency table.
#'
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param c0 numeric and greater than 0. It is a hyperparameter in "efron"
#' model.
#' @param pDegree integer greater than or equal to 2. It is a hyperparameter in
#' Efron model.
#' @param aStart initial value for parameter alpha in Efron model.
#'
#' @details
#'
#' The implementation of the "efron" model is adapted from the
#' \pkg{deconvolveR} package, developed by Bradley Efron and
#' Balasubramanian Narasimhan. The original implementation in \pkg{deconvolveR}
#' does not support an exposure or offset parameter in the Poisson model,
#' which corresponds to the expected null value (\eqn{E_{ij}}) for each AE-drug
#' combination. To address this, we modified the relevant code to allow for the
#' inclusion of \eqn{E_{ij}} in the Poisson likelihood. In addition, we
#' implemented a method for estimating the degrees of freedom, enabling AIC- or
#' BIC-based hyperparameter selection for the "efron" model (Tan et al. 2025).
#' See \code{\link{pvEBayes_tune}} for details.
#'
#' @references
#'
#' Narasimhan B, Efron B. deconvolveR: A G-modeling program for deconvolution
#' and empirical Bayes estimation. \emph{Journal of Statistical Software}.
#' 2020; 2;94:1-20.
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
#'
#'
#' @returns a list of optimizer outputs
#' @keywords internal
.E_fit <- function(N, E, c0 = 1, pDegree = 5, aStart = 1, rel.tol) {
  tau <- .grid_based_on_hist_log_scale_sobol(N, E, max_draws = 3010)
  x <- as.vector(N)
  E <- as.vector(E)
  P <- vapply(tau, function(ee) {
    stats::dpois(x, ee * E)
  }, FUN.VALUE = numeric(length(x)))
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
    qw <- crossprod(Q, W)
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
    qw <- crossprod(Q, W)
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
    hessian = hess_fun, control = list(rel.tol = rel.tol)
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
    qw <- crossprod(Q, W)
    qw_rowsums <- rowSums(qw)
    qg <- crossprod(Q, g)
    t1 <- tcrossprod(qw)
    t2 <- tcrossprod(qw_rowsums, qg)
    t3 <- t(t2)
    W_rowsums <- rowSums(W)
    t4 <- crossprod(Q * W_rowsums, Q)
    hess <- (t1 + t2 + t3 - t4)
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





#' Fit gamma mixture based empirical Bayes models using ECM algorithm.
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param dirichlet logical. Used for "general-gamma" model. If is TRUE, a
#' dirichlet hyperprior for weights of gamma mixture prior is applied.
#' @param alpha numeric between 0 and 1. The hyperparameter of "general-gamma"
#' model. It is needed if "general-gamma" model is used.
#' @param K integer greater than or equal to 2. It is needed if "K-gamma" model
#' is used.
#' @param maxi upper limit of iteration for the ECM algorithm.
#' @param h a vector of initialization of parameter h.
#' @param eps a tolerance parameter for ECM algorithm.
#'
#' @details
#' This function implements the ECM algorithm proposed by Tan et al. (2025),
#' providing a stable and efficient implementation of Gamma-Poisson
#' Shrinker(GPS), K-gamma and "general-gamma" methods for signal estimation and
#' signal detection in Spontaneous Reporting System (SRS) data table.
#'
#' @references
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
#'
#' DuMouchel W. Bayesian data mining in large frequency tables, with an
#' application to the FDA spontaneous reporting system.
#' \emph{The American Statistician.} 1999; 1;53(3):177-90. \cr
#'
#' @returns a list of optimizer outputs
#'
#' @keywords internal
#'
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
      message(
        paste0(
          "parameter alpha is not needed in GPS/K-gamma model,",
          " now running without alpha."
        )
      )
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
      message(
        paste0(
          "parameter K is not needed in general-gamma model,",
          " now running without K."
        )
      )
    }
    n_entry <- prod(dim(N))
    grid <- .grid_based_on_hist_log_scale_sobol(N, E,
      max_draws = min(200, n_entry)
    )
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
  result
}


#' Calculate mixture weights in posterior distribution for gamma mixture based
#' empirical Bayes models
#'
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param log logical. If is TRUE, log-sum-exp trick is used in the computation.
#'
#' @returns an IxJ by K matrix, where K is the number of mixture components.
#' @keywords internal
#' @noRd
.calc_qn_mix_gamma <- function(object, N, E, log = FALSE) {
  num_comp <- length(object$r)
  if (num_comp == 1) {
    res <- matrix(1, length(N), 1)
    return(res)
  }
  alphas <- object$r
  h <- object$h
  P <- object$omega
  N <- as.vector(N)
  E <- as.vector(E)
  if (log == FALSE) {
    Tau <- vapply(seq_along(N), function(e) {
      tmp <- stats::dnbinom(N[e],
        size = alphas,
        prob = 1 / (1 + E[e] * h), log = TRUE
      ) + log(P)
      tmp <- exp(tmp - max(tmp))
    }, FUN.VALUE = numeric(length(alphas))) %>% t()
    Qn <- Tau / rowSums(Tau)
  } else {
    Tau <- vapply(seq_along(N), function(e) {
      tmp <- stats::dnbinom(N[e],
        size = alphas,
        prob = 1 / (1 + E[e] * h), log = TRUE
      ) + log(P)
      tmp <- tmp - max(tmp)
    }, FUN.VALUE = numeric(length(alphas))) %>% t()
    Qn <- Tau - log(rowSums(exp(Tau)))
  }

  Qn
}

#' Take posterior draws for gamma mixture based models
#'
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param nsim number of posterior draws for each AE-drug combination.
#'
#' @returns a \code{pvEBayes} object
#' @keywords internal
#' @noRd
.generate_posterior_gamma_mix <- function(N,
                                          E,
                                          object,
                                          nsim = 1000) {
  prior_shape <- object$r
  prior_rate <- 1 / object$h
  m <- length(prior_shape)
  I <- nrow(N)
  J <- ncol(N)

  Qn <- .calc_qn_mix_gamma(object, N, E)
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

  lambda_sim
}


#' Calculate mixture weights in posterior distribution for grid based
#' empirical Bayes models
#'
#' @param N an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table.
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param nsim number of posterior draws for each AE-drug combination.
#'
#' @returns a \code{pvEBayes} object
#' @keywords internal
#' @noRd
.generate_posterior_grid_based <- function(N, E,
                                           object,
                                           nsim = 1000) {
  grid <- object$grid
  esti_prior <- object$g
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

  lambda_sim
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
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param model the model to fit. Available models are "general-gamma",
#' "K-gamma", "GPS", "KM" and "efron". Default to "general-gamma". Note that the
#' input for model is case-sensitive.
#' @param alpha numeric between 0 and 1. The hyperparameter of "general-gamma"
#' model. It is needed if "general-gamma" model is used. Small 'alpha'
#' encourages shrinkage on mixture weights of the estimated prior distribution.
#' See Tan et al. (2025) for further details.
#' @param K a integer greater than or equal to 2 indicating the number of
#' mixture components in the prior distribution. It is needed if "K-gamma"
#' model is used. See Tan et al. (2025) for further details.
#' @param p a integer greater than or equal to 2. It is needed if "efron" mode
#' is used. Larger p leads to smoother estimated prior distribution. See
#' Narasimhan and Efron (2020) for detail.
#' @param c0 numeric and greater than 0. It is needed if "efron" mode is used.
#' Large c0 encourage estimated prior distribution shrink toward discrete
#' uniform. See Narasimhan and Efron (2020) for detail.
#' @param maxi a upper limit of iteration for the ECM algorithm.
#' @param n_posterior_draws a number of posterior draws for each AE-drug
#' combination.
#' @param tol_ecm a tolerance parameter used for the ECM stopping rule, defined
#' as the absolute change in the joint marginal likelihood between two
#' consecutive iterations. It is used when 'GPS', 'K-gamma' or 'general-gamma'
#' model is fitted. Default to be 1e-4.
#' @param rtol_efron a tolerance parameter used when 'efron' model is fitted.
#' Default to 1e-10. See 'stats::nlminb' for detail.
#' @param rtol_KM a tolerance parameter used when 'KM' model is fitted.
#' Default to be 1e-6. See 'CVXR::solve' for detail.
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table. If `NULL` (default), the expected counts are estimated
#' from the SRS data using 'estimate_null_expected_count()'.
#'
#'
#' @details
#'
#' This function implements the ECM algorithm proposed by Tan et al. (2025),
#' providing a stable and efficient implementation of Gamma-Poisson
#' Shrinker(GPS), K-gamma and "general-gamma" methods for signal estimation and
#' signal detection in Spontaneous Reporting System (SRS) data table.
#'
#' Method "GPS" is proposed by DuMouchel (1999) and it is a parametric empirical
#' Bayes model with a two gamma mixture prior distribution.
#'
#' Methods "K-gamma" and "general-gamma" are non-parametric empirical Bayes
#' models, introduced by Tan et al. (2025). The number of mixture components "K"
#' needs to be prespecified when fitting a "K-gamma" model. For "general-gamma",
#' the mixture weights are regularized by a Dirichlet hyper prior with
#' hyperparameter \eqn{0 \leq \alpha < 1} that controls the shrinkage strength.
#' As "alpha" approaches 0, less non-empty mixture components exist in the
#' fitted model. When \eqn{\alpha = 0}, the Dirichlet distribution is an
#' improper prior still offering a reasonable posterior inference that
#' represents the strongest shrinkage of the "general-gamma" model.
#'
#' Parameter estimation for the "KM" model is formulated as a convex
#' optimization problem. The objective function and constraints follow the same
#' construction as in \pkg{REBayes}. Parameter estimation is performed using
#' the open-source convex optimization package \pkg{CVXR}.
#'
#'
#' The implementation of the "efron" model in this package is adapted from the
#' \pkg{deconvolveR} package, developed by Bradley Efron and
#' Balasubramanian Narasimhan. The original implementation in \pkg{deconvolveR}
#' does not support an exposure or offset parameter in the Poisson model,
#' which corresponds to the expected null value (\eqn{E_{ij}}) for each AE-drug
#' combination. To address this, we modified the relevant code to allow for the
#' inclusion of \eqn{E_{ij}} in the Poisson likelihood. In addition, we
#' implemented a method for estimating the degrees of freedom, enabling AIC- or
#' BIC-based hyperparameter selection for the "efron" model (Tan et al. 2025).
#' See \code{\link{pvEBayes_tune}} for details.
#'
#'
#'
#' @references
#'
#' DuMouchel W. Bayesian data mining in large frequency tables, with an
#' application to the FDA spontaneous reporting system.
#' \emph{The American Statistician.} 1999; 1;53(3):177-90. \cr
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
#'
#' Narasimhan B, Efron B. deconvolveR: A G-modeling program for deconvolution
#' and empirical Bayes estimation. \emph{Journal of Statistical Software}.
#' 2020; 2;94:1-20.
#'
#' Koenker R, Gu J. REBayes: an R package for empirical Bayes mixture methods.
#' \emph{Journal of Statistical Software}. 2017; 4;82:1-26.
#'
#' Fu, A, Narasimhan, B, Boyd, S. CVXR: An R Package for Disciplined Convex
#' Optimization. \emph{Journal of Statistical Software}. 2020; 94;14:1-34.
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
#' sig_mat[c(1, 5), 1] <- 2
#'
#'
#' simu_table <- generate_contin_table(
#'   ref_table = ref_table,
#'   signal_mat = sig_mat,
#'   n_table = 1
#' )[[1]]
#'
#'
#' # fit general-gamma model with a specified alpha
#' fit <- pvEBayes(
#'   contin_table = simu_table,
#'   model = "general-gamma",
#'   alpha = 0.3,
#'   n_posterior_draws = 1000,
#'   E = NULL,
#'   maxi = NULL
#' )
#'
#' # fit K-gamma model with K = 3
#' fit_Kgamma <- pvEBayes(
#'   contin_table = simu_table, model = "K-gamma",
#'   K = 3, n_posterior_draws = 1000
#' )
#'
#'
#' # fit Efron model with specified hyperparameters
#' # p = 40, c0 = 0.05
#'
#' \dontrun{
#' fit_efron <- pvEBayes(
#'   contin_table = simu_table,
#'   model = "efron",
#'   p = 40,
#'   c0 = 0.05,
#'   n_posterior_draws = 1000
#' )
#' }
#'
#' # fit GPS model and comapre with 'openEBGM'
#'
#'
#' fit_gps <- pvEBayes(simu_table, model = "GPS")
#'
#' \dontrun{
#'
#' ## Optional comparison with openEBGM (only if installed)
#'
#' ## tol_ecm is the absolute tolerance for ECM stopping rule.
#' ## It is set to ensure comparability to `openEBGM`.
#'
#' fit_gps <- pvEBayes(simu_table, model = "GPS", tol_ecm = 1e-2)
#'
#' if (requireNamespace("openEBGM", quietly = TRUE)) {
#'   E <- estimate_null_expected_count(simu_table)
#'   simu_table_stacked <- as.data.frame(as.table(simu_table))
#'   simu_table_stacked$E <- as.vector(E)
#'   colnames(simu_table_stacked) <- c("var1", "var2", "N", "E")
#'   simu_table_stacked_squash <- openEBGM::autoSquash(simu_table_stacked)
#'
#'   hyper_estimates <- openEBGM::hyperEM(simu_table_stacked_squash,
#'     theta_init = c(2, 1, 2, 2, 0.2),
#'     method = "nlminb",
#'     N_star = NULL,
#'     zeroes = TRUE,
#'     param_upper = Inf,
#'     LL_tol = 1e-2,
#'     max_iter = 10000
#'   )
#' }
#'
#' theta_hat <- hyper_estimates$estimates
#' qn <- openEBGM::Qn(theta_hat,
#'   N = simu_table_stacked$N,
#'   E = simu_table_stacked$E
#' )
#'
#' simu_table_stacked$q05 <- openEBGM::quantBisect(5,
#'   theta_hat = theta_hat,
#'   N = simu_table_stacked$N,
#'   E = simu_table_stacked$E,
#'   qn = qn
#' )
#'
#' ## obtain the detected signal provided by openEBGM
#' simu_table_stacked %>%
#'   subset(q05 > 1.001)
#'
#' ## detected signal from pvEBayes presented in the same way as openEBGM
#' fit_gps %>%
#'   summary(return = "posterior draws") %>%
#'   apply(c(2, 3), quantile, prob = 0.05) %>%
#'   as.table() %>%
#'   as.data.frame() %>%
#'   subset(Freq > 1.001)
#' }
#'
#' @srrstats {G1.6} We provide an implementation of GPS model with the ECM
#' algorithm as an alternative to \code{openEBGM} package. A comparison of
#' these two implementation is presented in examples of associated function
#' documentation.
#' @srrstats {G2.0, G2.1, G2.2} Length and value of single and vector inputs are
#' properly checked.
#' @srrstats {G2.0a, G2.1a} The length of single and vector inputs are
#' explicitly described in the corresponding documentation.
#' @srrstats {G2.3, G2.3a, G2.3b} Character inputs are explicitly documented
#' that they are strictly case-sensitive and only applicable to expected values.
#' @srrstats {G2.4, G2.4a, G2.4b, G2.4c, G2.8} Explicit conversion is used for
#' integer, continuous and character inputs.
#' @srrstats {G2.7} Tabular formats appear in Depends or Sugggests are tested.
#' @srrstats {G3.0} The algorithm do not compare floating points for equality.
#' @srrstats {G5.4b} The examples field provide a comparison of implementation
#' of GPS model showing that both implementation correctly detect the signal.
#' @srrstats {BS5.3, BS5.4, BS5.5}
#' The empirical Bayes methods implemented in \pkg{pvEBayes} do not rely on
#' stochastic sampling, and therefore do not produce the types of
#' convergence diagnostics typically associated with full Bayesian modeling.
#' Convergence in the ECM algorithm is reached (at least to a sub-optimal).
#' This is ensured by monotonically increased log joint marginal likelihood,
#' as proved by Tan et al. (*Stat. in Med*, 2025).
#'
pvEBayes <- function(contin_table, model = "general-gamma",
                     alpha = NULL, K = NULL,
                     p = NULL, c0 = NULL,
                     maxi = NULL,
                     tol_ecm = 1e-4,
                     rtol_efron = 1e-10,
                     rtol_KM = 1e-6,
                     n_posterior_draws = 1000,
                     E = NULL) {
  h <- NULL
  contin_table <- as.matrix(contin_table)
  if (.is_valid_contin_table(contin_table) == FALSE) {
    stop()
  }
  if (is.null(colnames(contin_table)) ||
    is.null(rownames(contin_table))) {
    contin_table <- .set_default_names(contin_table)
  }
  if (is.null(E)) {
    E <- calculate_tilde_e(contin_table)
  } else {
    if (!(all(E >= 0) &&
      identical(dim(E), dim(contin_table)))
    ) {
      stop(
        paste0(
          "'E' must contain only positive values and have the same ",
          "dimensions as 'contin_table'."
        )
      )
    }
  }
  E <- as.matrix(E)
  if (!(length(model) == 1 &&
    (model %in% c("general-gamma", "K-gamma", "GPS", "KM", "efron")))) {
    stop(
      paste0(
        "'model' must be one of the followings:",
        "'general-gamma', 'K-gamma', 'GPS', 'KM', 'efron'"
      )
    )
  }
  model <- as.character(model)
  if (!is.null(maxi) &&
    !(is.numeric(maxi) && length(maxi) == 1 &&
      maxi %% 1 == 0 && maxi > 0)) {
    stop("'maxi' must be a single integer that is greater than 0.")
  }
  if (!is.null(maxi)) {
    maxi <- as.integer(maxi)
  } else {
    maxi <- 1000L
  }

  if (!is.null(n_posterior_draws) &&
    !(is.numeric(n_posterior_draws) && length(n_posterior_draws) == 1 &&
      n_posterior_draws %% 1 == 0 && n_posterior_draws > 0)) {
    stop("'n_posterior_draws' must be a single positive integer.")
  }

  if (!is.null(n_posterior_draws)) {
    n_posterior_draws <- as.integer(n_posterior_draws)
  }


  if (!(is.numeric(tol_ecm) && length(tol_ecm) == 1 &&
    tol_ecm > 0)) {
    stop("'tol_ecm' must be a single positive variable.")
  }
  tol_ecm <- as.numeric(tol_ecm)

  if (model == "general-gamma") {
    if (!(is.numeric(alpha) && length(alpha) == 1 &&
      alpha >= 0 && alpha <= 1)) {
      stop("'alpha' must be a single numeric variable between 0 and 1.")
    }
    alpha <- as.numeric(alpha)
    dirichlet <- TRUE
  } else {
    if (model == "GPS") {
      K <- 2
      dirichlet <- FALSE
    }
    if (model == "K-gamma") {
      if (!(is.numeric(K) && length(K) == 1 &&
        K %% 1 == 0 && K > 2)) {
        stop("'K' must be a single integer that is greater than 2.")
      }
      K <- as.integer(K)
      dirichlet <- FALSE
    }
  }
  start_time <- Sys.time()
  if (model == "KM") {
    res <- .KM_fit(
      contin_table,
      E,
      rtol_KM = rtol_KM
    )
    generate_posterior_fun <- .generate_posterior_grid_based
  } else if (model == "efron") {
    if (!(is.numeric(p) && length(p) == 1 &&
      p %% 1 == 0 && p >= 2)) {
      stop("'p' must be a single integer that is greater or equal to 2.")
    }
    p <- as.integer(p)

    if (!(is.numeric(c0) && length(c0) == 1 &&
      c0 > 0)) {
      stop("'c0' must be a single numberic variable that is greater than 0.")
    }
    c0 <- as.numeric(c0)
    res <- .E_fit(
      contin_table,
      E,
      c0,
      p,
      rel.tol = rtol_efron
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
      tol_ecm
    )
    res$alpha <- alpha
    generate_posterior_fun <- .generate_posterior_gamma_mix
  }

  end_time <- Sys.time() # start the clock

  res$fit_time <- difftime(end_time, start_time)
  res$draws_time <- NULL
  res$model <- model


  if (is.null(n_posterior_draws)) {
    res$n_posterior_draws <- n_posterior_draws
    res$draws_time <- NULL
  } else {
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
#' function fits a series of models across a grid of candidate hyperparameter
#' values and computes their AIC and BIC. The models with the lowest AIC or BIC
#' values are selected as the optimal fits.
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param model the model to be tuned. Available models are "general-gamma" and
#' "efron". Default to "general-gamma". Note that the input for model is
#' case-sensitive.
#' @param alpha_vec vector of hyperparameter alpha values to be selected. Alpha
#' is a hyperparameter in general-gamma model which is numeric and between 0
#' and 1. If is NULL, a default set of alpha values (0, 0.1, 0.3, 0.5, 0.7, 0.9)
#' will be used.
#' @param p_vec vector of hyperparameter p values to be selected. p is
#' a hyperparameter in "efron" model which should be a positive integer. If is
#' NULL, a default set of p values (80, 100, 120, 150, 200) will be used.
#' @param c0_vec vector of hyperparameter c0 values to be selected. c0 is
#' a hyperparameter in "efron" model which should be a positive number. If is
#' NULL, a default set of c0 values (0.001, 0.01, 0.1, 0.2, 0.5) will be used.
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
#' @param tol_ecm a tolerance parameter used for the ECM stopping rule, defined
#' as the absolute change in the joint marginal likelihood between two
#' consecutive iterations. It is used when 'GPS', 'K-gamma' or 'general-gamma'
#' model is fitted. Default to be 1e-4.
#' @param rtol_efron a tolerance parameter used when 'efron' model is fitted.
#' Default to 1e-10. See 'stats::nlminb' for detail.
#' @param E A matrix of expected counts under the null model for the SRS
#' frequency table. If `NULL` (default), the expected counts are estimated
#' from the SRS data using 'estimate_null_expected_count()'.
#'
#'
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control.} 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model.
#' \emph{The Annals of Statistics.} 1978; 1:461-4.
#'
#' Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
#' to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
#' Estimation in Spontaneous Reporting Systems Data.
#' \emph{Statistics in Medicine.} 2025; 44: 18-19,
#' https://doi.org/10.1002/sim.70195.
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
#' fit <- pvEBayes_tune(statin2025_44,
#'   model = "general-gamma",
#'   alpha_vec = c(0, 0.1, 0.3, 0.5, 0.7, 0.9)
#' )
#'
#' @srrstats {G2.0, G2.1, G2.2} length and value of single and
#' vector inputs are properly checked.
#' @srrstats {G2.0a, G2.1a} The length of single and vector inputs are explicitly
#' described in the corresponding documentation.
#' @srrstats {G2.3, G2.3a, G2.3b} character inputs are explicitly documented
#' that they are strictly case-sensitive and only applicable to expected values.
#' @srrstats {G2.4, G2.4a, G2.4b, G2.4c, G2.8} explicit conversion is used for
#' integer, continuous and character inputs.
#' @srrstats {G2.6} one-dimensinoal inputs are appropriately pre-processed.
pvEBayes_tune <- function(contin_table, model = "general-gamma",
                          alpha_vec = NULL,
                          p_vec = NULL, c0_vec = NULL,
                          use_AIC = TRUE,
                          n_posterior_draws = 1000,
                          return_all_fit = FALSE,
                          return_all_AIC = TRUE,
                          return_all_BIC = TRUE,
                          tol_ecm = 1e-4,
                          rtol_efron = 1e-10,
                          E = NULL) {
  if (!.is_valid_contin_table(contin_table)) {
    stop("Please provide a valid 'contin_table'")
  }
  if (!(is.logical(return_all_fit) && length(return_all_fit) == 1)) {
    stop("'return_all_fit' must be a single logical value (TRUE or FALSE).")
  }

  if (!(is.logical(return_all_AIC) && length(return_all_AIC) == 1)) {
    stop("'return_all_AIC' must be a single logical value (TRUE or FALSE).")
  }

  if (!(is.logical(return_all_BIC) && length(return_all_BIC) == 1)) {
    stop("'return_all_BIC' must be a single logical value (TRUE or FALSE).")
  }

  if (is.null(colnames(contin_table)) ||
    is.null(rownames(contin_table))) {
    contin_table <- .set_default_names(contin_table)
  }
  if (is.null(E)) {
    E <- calculate_tilde_e(contin_table)
  } else {
    if (!(all(E >= 0) &&
      identical(dim(E), dim(contin_table)))
    ) {
      stop(
        paste0(
          "'E' must contain only positive values and have the same ",
          "dimensions as 'contin_table'."
        )
      )
    }
  }
  E <- as.matrix(E)
  if (!(length(model) == 1 &&
    (model %in% c("general-gamma", "K-gamma", "GPS", "KM", "efron")))) {
    stop(
      paste0(
        "'model' must be one of the followings:",
        "'general-gamma', 'K-gamma', 'GPS', 'KM', 'efron'"
      )
    )
  }
  model <- as.character(model)
  if (model %in% c("K-gamma", "GPS", "KM")) {
    stop("Please use pvEBayes() for GPS, K-gamma or KM model fitting.")
  }

  if (!is.null(n_posterior_draws) &&
    !(is.numeric(n_posterior_draws) && length(n_posterior_draws) == 1 &&
      n_posterior_draws %% 1 == 0 && n_posterior_draws > 0)) {
    stop("'n_posterior_draws' must be a single positive integer.")
  }

  if (model == "efron") {
    if (is.null(p_vec)) {
      p_vec <- c(40, 60, 80, 100, 120)
    }
    if (is.null(c0_vec)) {
      c0_vec <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
    }
    p_vec <- as.vector(p_vec)
    c0_vec <- as.vector(c0_vec)
    if (!(
      length(p_vec) > 0 &&
        length(c0_vec) > 0 &&
        all(p_vec >= 2) &&
        all(c0_vec >= 0)
    )
    ) {
      stop("Please provide valid 'p_vec' and 'c0_vec'.")
    }
    p_vec <- as.integer(p_vec)
    c0_vec <- as.numeric(c0_vec)
    objects <- tuning_efron(
      contin_table = contin_table,
      p_vec = p_vec,
      c0_vec = c0_vec,
      return_all_fit = return_all_fit,
      return_all_AIC = return_all_AIC,
      return_all_BIC = return_all_BIC,
      rtol_efron = rtol_efron
    )
    generate_posterior_fun <- .generate_posterior_grid_based
  } else {
    if (is.null(alpha_vec)) {
      alpha_vec <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9)
    }
    alpha_vec <- as.vector(alpha_vec)
    if (!(
      length(alpha_vec) > 0 &&
        is.numeric(alpha_vec) &&
        all(alpha_vec >= 0 & alpha_vec <= 1)
    )
    ) {
      stop("Elements in 'alpha_vec' must be numeric and in [0,1)")
    }
    alpha_vec <- as.numeric(alpha_vec)
    objects <- tuning_general_gamma(
      contin_table = contin_table,
      alpha_vec = alpha_vec,
      return_all_fit = return_all_fit,
      return_all_AIC = return_all_AIC,
      return_all_BIC = return_all_BIC,
      tol_ecm = tol_ecm
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


  res$posterior_draws <- generate_posterior_fun(contin_table,
    E,
    res,
    nsim = n_posterior_draws
  )
  end_time2 <- Sys.time()
  res$draws_time <- difftime(end_time2, end_time)

  res$contin_table <- contin_table
  res$E <- E

  res$tuning <- objects

  class(res) <- c("pvEBayes", "pvEBayes_tuned")
  res
}




#' Select hyperparameter alpha and obtain the optimal general-gamma model based
#'  on AIC and BIC
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param alpha_vec vector of hyperparameter alpha values to be selected. Alpha
#' is hyperparameter in general-gamma model which is numeric and between 0
#' and 1. If is NULL, a default set of alpha values (0, 0.1, 0.3, 0.5, 0.7, 0.9)
#' will be used.
#' @param tol_ecm a tolerance parameter used for the ECM stopping rule, defined
#' as the absolute change in the joint marginal likelihood between two
#' consecutive iterations. It is used when 'GPS', 'K-gamma' or 'general-gamma'
#' model is fitted. Default to be 1e-4.
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control.} 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model.
#' \emph{The Annals of Statistics.} 1978; 1:461-4.
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
                                 return_all_BIC = TRUE,
                                 tol_ecm = 1e-4) {
  fits <- alpha_vec %>%
    lapply(function(e) {
      pvEBayes(
        contin_table = contin_table, model = "general-gamma",
        alpha = e, n_posterior_draws = NULL,
        tol_ecm = tol_ecm
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


#' generate a list of all possible hyperparameter combination for Efron model
#'
#' @param p a vector of hyperparameter p values to be selected.
#' @param c0 a vector of hyperparameter c0 values to be selected.
#'
#' @returns a list
#' @keywords internal
#' @noRd
.get_grid_list <- function(c0, p) {
  res <- list()
  for (i in seq_along(c0)) {
    for (j in seq_along(p)) {
      res[[10 * i + j]] <- list(c0 = c0[i], p = p[j])
    }
  }
  res <- Filter(Negate(is.null), res)
  res
}


#' Select hyperparameter (p, c0) and obtain the optimal efron model based
#'  on AIC and BIC
#'
#'
#' @param contin_table an IxJ contingency table showing pairwise counts of
#' adverse events for I AEs (along the rows) and J drugs (along the columns).
#' @param p_vec vector of hyperparameter p values to be selected. p is
#' a hyperparameter in "efron" model which should be a positive integer. If is
#' NULL, a default set of p values (80, 100, 120, 150, 200) will be used.
#' @param c0_vec vector of hyperparameter c0 values to be selected. c0 is
#' a hyperparameter in "efron" model which should be a positive number. If is
#' NULL, a default set of c0 values (0.001, 0.01, 0.1, 0.2, 0.5) will be used.
#' @param rtol_efron a tolerance parameter used when 'efron' model is fitted.
#' Default to 1e-10. See 'stats::nlminb' for detail.
#'
#' @references
#'
#' Akaike H. A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control.}
#' 2003; 19(6):716-23. \cr
#'
#' Schwarz G. Estimating the dimension of a model.
#' \emph{The Annals of Statistics.} 1978; 1:461-4.
#'
#' @return
#' a list of fitted models with hyperparameter alpha selected by AIC or BIC.
#'
#' @keywords internal
tuning_efron <- function(contin_table,
                         p_vec = NULL,
                         c0_vec = NULL,
                         return_all_fit = FALSE,
                         return_all_AIC = TRUE,
                         return_all_BIC = TRUE,
                         rtol_efron = 1e-10) {
  parameter_list <- .get_grid_list(c0_vec, p_vec)
  fits <- parameter_list %>%
    lapply(function(e) {
      pvEBayes(
        contin_table = contin_table, model = "efron",
        p = e$p, c0 = e$c0,
        n_posterior_draws = NULL,
        rtol_efron = rtol_efron
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
