#' Generate random sample from a truncated normal distribution
#'
#' @param n number of draws.
#' @param a lower bound of the truncation.
#' @param b upper bound of the truncation.
#' @param mean mean of the normal distribution.
#' @param sd standard deviation of the normal distribution.
#'
#' @returns a vector of draws
#' @keywords internal
#' @noRd
rtruncnorm <- function(n, a = -0.4, b = 0, mean = 0, sd = 0.05) {
  samples <- list()
  count <- 0 # Track the number of accepted samples
  i <- 1
  while (count < n) {
    # Generate more than needed to reduce iterations
    m <- ceiling(1.2 * n)
    candidates <- stats::rnorm(m, mean = mean, sd = sd)

    # Keep only the values within the truncation bounds
    samples[[i]] <- candidates[(candidates >= a) & (candidates <= b)]
    count <- count + length(samples[[i]])
    # Store in preallocated vector
    i <- i + 1
  }
  res <- unlist(samples)[1:n]
  res
}




#' Generate random SRS tables
#'
#' @param ref_table a reference table used as the basis for generating random
#' tables.
#' @param signal_mat numeric matrix of the same dimension as
#' the reference table (ref_table). The entry at position (i, j)
#' in signal_mat represents the signal strength between
#' the i-th adverse event and the j-th drug. By default,
#' each pair is assigned a value of 1, indicating no signal for that pair.
#' @param n_simu number of random matrices to generate.
#' @param Variation logical. Include random noises to sig_mat while
#' generating random tables. Default to FALSE.
#' If set to TRUE, n_table of sig_mat incorporating the added noise
#' will also be returned.
#' @param Z logical matrix of the same size as ref_table indicating
#' the positions of structural zero.
#'
#' @returns
#'
#' A list of length \code{n_table}, with each entry being a
#' \code{nrow(ref_table)} by \code{ncol(ref_table)} matrix.
#'
#' @keywords internal
#' @noRd
.simu_table_multinom <- function(ref_table,
                                 signal_mat,
                                 n_simu,
                                 Variation = FALSE,
                                 Z = NULL) {
  Ndd <- sum(ref_table)
  pid <- rowSums(ref_table) / Ndd
  pdj <- colSums(ref_table) / Ndd
  I <- nrow(ref_table)
  J <- ncol(ref_table)
  if (is.null(Z)) {
    Z <- matrix(0, I, J)
  }
  a <- sum(pid %*% t(pdj) * signal_mat * (1 - Z))
  if (Variation == TRUE) {
    var_sig <- lapply(1:n_simu, function(e) {
      tmp_sig <- signal_mat
      err1 <- matrix(rtruncnorm(
        I * J,
        a = -0.4,
        b = 0,
        mean = 0,
        sd = 0.05
      ), I, J)
      err2 <- matrix(rtruncnorm(
        I * J,
        a = -0.2,
        b = 0.2,
        mean = 0,
        sd = 0.05
      ), I, J)
      tmp_sig[signal_mat == 1] <- signal_mat[signal_mat == 1] +
        err1[signal_mat == 1]
      tmp_sig[signal_mat > 1] <- signal_mat[signal_mat > 1] +
        err2[signal_mat > 1]
      tmp_sig[I, ] <- 1
      tmp_sig[, J] <- 1
      tmp_sig * (1 - Z) %>%
        magrittr::set_rownames(rownames(ref_table)) %>%
        magrittr::set_colnames(colnames(ref_table))
    })
  }
  tables <- lapply(1:n_simu, function(e) {
    if (Variation == TRUE) {
      tmp_sig <- var_sig[[e]]
      a <- sum(pid %*% t(pdj) * tmp_sig)
    } else {
      tmp_sig <- signal_mat * (1 - Z)
    }
    matrix(stats::rmultinom(
      n = 1,
      size = Ndd,
      prob = pid %*% t(pdj) * tmp_sig / a
    ), I, J) %>%
      magrittr::set_rownames(rownames(ref_table)) %>%
      magrittr::set_colnames(colnames(ref_table))
  })
  if (Variation == TRUE) {
    list(
      tables = tables,
      var_sig = var_sig
    )
  } else {
    tables
  }
}


#' Generate random contingency tables based on a reference table
#' embedded signals,and possibly with zero inflation
#'
#' @description
#' This function generates random contingency tables that resemble a given
#' reference table, with the option to embed signals and zero-inflation.
#'
#'
#' @param n_table a number of random matrices to generate.
#' @param ref_table a reference table used as the basis for generating random
#' tables.
#' @param signal_mat a numeric matrix of the same dimension as
#' the reference table (ref_table). The entry at position (i, j)
#' in signal_mat represents the signal strength between
#' the i-th adverse event and the j-th drug. By default,
#' each pair is assigned a value of 1, indicating no signal for that pair.
#' @param Variation logical. Include random noises to sig_mat while
#' generating random tables. Default to FALSE.
#' If set to TRUE, n_table of sig_mat incorporating the added noise
#' will also be returned.
#' @param zi_indic_mat logical matrix of the same size as ref_table indicating
#' the positions of structural zero.
#'
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
#' A list of length \code{n_table}, with each entry being a
#' \code{nrow(ref_table)} by \code{ncol(ref_table)} matrix.
#'
#' @examples
#'
#' set.seed(1)
#' ref_table <- statin2025_44
#'
#'
#' # set up signal matrix with one signal at entry (1,1)
#' sig_mat <- matrix(1, nrow(ref_table), ncol(ref_table))
#' sig_mat[1, 1] <- 2
#'
#' # set up structural zero matrix
#' Z <- matrix(0, nrow(ref_table), ncol(ref_table))
#' Z[5, 1] <- 1
#'
#' simu_table <- generate_contin_table(ref_table,
#'   signal_mat = sig_mat,
#'   n_table = 1,
#'   Variation = TRUE,
#'   zi_indic_mat = Z
#' )[[1]][[1]]
#'
#' @export
#'
#' @srrstats {G2.0, G2.1, G2.2} length and value of "n_table" and "Variation"
#' are properly checked.
#' @srrstats {G2.0a, G2.1a} The length of "n_table" and "Variation" are explicitly
#' described in the corresponding documentation.
#' @srrstats {G2.4, G2.4a, G2.8} explicit conversion is used for integer input.
#'
generate_contin_table <- function(n_table = 1,
                                  ref_table,
                                  signal_mat = NULL,
                                  Variation = FALSE,
                                  zi_indic_mat = NULL) {
  if (is.null(signal_mat)) {
    signal_mat <- ref_table
    signal_mat[] <- 1
  }


  if (!(is.numeric(n_table) && length(n_table) == 1 &&
    n_table %% 1 == 0 && n_table > 0)) {
    stop("'n_table' must be a single integer.")
  }
  n_table <- as.integer(n_table)
  if (!(is.logical(Variation) && length(Variation) == 1)) {
    stop("'Variation' must be a single logical value (TRUE or FALSE).")
  }

  stopifnot(
    is.numeric(n_table),
    n_table > 0,
    nrow(ref_table) == nrow(signal_mat),
    ncol(ref_table) == ncol(signal_mat)
  )
  if (!is.null(zi_indic_mat)) {
    stopifnot(
      nrow(ref_table) == nrow(zi_indic_mat),
      ncol(ref_table) == ncol(zi_indic_mat)
    )
  } else {
    zi_indic_mat <- ref_table
    zi_indic_mat[] <- 0
  }



  out <- .simu_table_multinom(
    ref_table,
    signal_mat,
    n_table,
    Variation,
    zi_indic_mat
  )
  out
}
