is.pvEBayes <- function(object) {
  methods::is(object, "pvEBayes")
}

is.pvEBayes_tuned <- function(object) {
  methods::is(object, "pvEBayes_tuned")
}


#' Extract all fitted models from a tuned pvEBayes Object
#'
#' @description
#' This function retrieves the list of all fitted models from a pvEBayes_tuned
#' object, which is the output of the `pvEBayes_tune()` function.
#'
#'
#' @param object An object of class \code{pvEBayes_tuned}, usually returned by
#' \code{\link{pvEBayes_tune}}. This function will throw an error if the input
#' is not of the correct class.
#'
#' @returns
#' A list containing the results of each model fitted during the tuning process.
#'
#' @export
#'
#' @examples
#'
#' valid_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2)
#' rownames(valid_matrix) <- c("AE_1", "AE_2")
#' colnames(valid_matrix) <- c("drug_1", "drug_2", "drug_3", "drug_4")
#'
#' tuned_object = pvEBayes_tune(valid_matrix,
#'                             model = "general-gamma",
#'                             return_all_fit = TRUE)
#' extract_all_fitted_models(tuned_object)
#'
extract_all_fitted_models <- function(object) {
  if (!is.pvEBayes_tuned(object)) {
    stop("This function can only be used after tuning. Please apply to objects returned by 'pvEBayes_tune()'.")
  }
  return(object$tuning)
}




#' Generate posterior draws for each AE-drug combination
#'
#' @description
#' This function generates posterior draws from the posterior distribution of
#' \eqn{\lambda_{ij}} for each AE-drug combination, based on a fitted empirical
#' Bayes model. The posterior draws can be used to compute credible intervals,
#' visualize posterior distributions, or support downstream inference.
#'
#'
#' @param obj a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param n_posterior_draws number of posterior draws for each AE-drug
#' combination.
#'
#' @return
#'
#' The function returns an S3 object of class `pvEBayes` with posterior draws.
#'
#' @export
#'
#' @examples
#'
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = NULL
#' )
#'
#' fit_with_draws <- posterior_draws(fit, n_posterior_draws = 1000)
#'
posterior_draws <- function(obj,
                                     n_posterior_draws = 1000) {
  stopifnot(is.pvEBayes(obj))

  if (obj$model %in% c("KM", "efron")) {
    generate_posterior_fun <- .generate_posterior_grid_based
  } else {
    generate_posterior_fun <- .generate_posterior_gamma_mix
  }
  start_time <- Sys.time()
  obj$posterior_draws <- generate_posterior_fun(obj$contin_table,
    obj$E, obj,
    nsim = n_posterior_draws
  )
  end_time <- Sys.time()
  obj$fit_time <- difftime(end_time, start_time)
  obj
}


.rownames_to_column <- function(df, var = "rowname") {
  # Create a new column from row names
  df[[var]] <- row.names(df)

  # Reset row names to sequential integers
  row.names(df) <- NULL

  # Reorder columns to place the new column first
  df <- df[, c(var, setdiff(names(df), var))]

  return(df)
}

.left_join_base <- function(x, y, by) {
  merge(x, y, by = by, all.x = TRUE)
}


.get_posterior_prob <- function(obj) {
  tmp <- obj$posterior_draws
  # posterior::draws_of()
  (tmp > 1.001) %>%
    apply(c(2, 3), mean)
}

.check_AEs <- function(post_matrix, N, N_threshold) {
  I <- nrow(N)
  J <- ncol(N)
  indi <- rep(FALSE, I)
  for (i in seq_len(I)) {
    tmp_indi <- (post_matrix[i, ] == max(post_matrix[i, ]))
    if (sum(tmp_indi) > 0) {
      if (sum(N[i, tmp_indi] > N_threshold) > 0) {
        indi[i] <- TRUE
      }
    }
  }
  indi
}



#' Generate an eyeplot showing the distribution of posterior draws for
#' selected drugs and adverse events
#'
#' @description
#' This function creates an eyeplot to visualize the posterior distributions of
#' \eqn{\lambda_{ij}} for selected AEs and drugs. The plot displays
#' posterior median, 90 percent credible interval for each selected AE-drug
#' combination.
#'
#'
#' @param x a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param num_top_AEs number of most significant AEs appearing in the plot. Default to
#' 10.
#' @param num_top_drugs number of most significant drugs appearing in the plot. Default
#' to 7.
#' @param specified_AEs a vector of AE names that are specified to appear in the plot. If
#' a vector of AEs is given, argument num_top_AEs will be ignored.
#' @param specified_drugs a vector of drug names that are specified to appear in the plot.
#' If a vector of drugs is given, argument num_top_drugs will be ignored.
#' @param N_threshold integer greater than 0. Any AE-drug combination with
#' observation smaller than N_threshold will be filtered out.
#' @param text_shift numeric. Controls the relative position of text labels,
#' (e.g., "N = 1", "E = 2"). A larger value shifts the "E = 2" further away from
#' its original position.
#' @param x_lim_scalar numeric. An x-axis range scalar that ensures text labels
#' are appropriately included in the plot.
#' @param text_size numeric. Controls the size of text labels,
#' (e.g., "N = 1", "E = 2").
#' @param log_scale logical. If TRUE, the eye plot displays the posterior
#' distribution of \eqn{\log(\lambda_{ij})} for the selected AEs and drugs.
#'
#' @return
#' a ggplot2 object.
#' @export
#'
#' @examples
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = 1000
#' )
#'
#' AE_names <- rownames(statin2025_44)[1:6]
#' drug_names <- colnames(statin2025_44)[-7]
#'
#' eyeplot_pvEBayes(
#'   x = fit
#' )
#'
eyeplot_pvEBayes <- function(x,
                             num_top_AEs = 10,
                             num_top_drugs = 8,
                             specified_AEs = NULL,
                             specified_drugs = NULL,
                             N_threshold = 1,
                             text_shift = 4,
                             x_lim_scalar = 1.3,
                             text_size = 3,
                             log_scale = FALSE) {
  top_drugs <- num_top_drugs
  top_AEs <- num_top_AEs
  AEs <- specified_AEs
  drugs <- specified_drugs
  if (top_drugs > (ncol(x$E) - 1)) {
    top_drugs <- (ncol(x$E) - 1)
  }
  stopifnot(is.pvEBayes(x))
  if (is.null(x$posterior_draws)) {
    x <- x %>% posterior_draws
  }

  counts_long <- x$contin_table %>%
    as.data.frame() %>%
    .rownames_to_column(var = "AE") %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars = "AE", variable.name = "drug", value.name = "N")
  Es_long <- x$E


  Es_long <- Es_long %>%
    round(2) %>%
    as.data.frame() %>%
    .rownames_to_column(var = "AE") %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars = "AE", variable.name = "drug", value.name = "E")

  counts_long <- counts_long %>%
    .left_join_base(Es_long, by = c("AE", "drug"))


  post_prob_matrix <- x$posterior_draws %>%
    # posterior::draws_of() %>%
    {
      . > 1.001
    } %>%
    apply(c(2, 3), mean)
  RMSE1 <- x$posterior_draws %>%
    # posterior::draws_of() %>%
    {
      (. - 1)^2
    } %>%
    apply(c(2, 3), mean)
  filter_indi <- .check_AEs(post_prob_matrix, x$contin_table, N_threshold)
  RMSE1[x$contin_table <= N_threshold] <- 0
  orders <- RMSE1[filter_indi, , drop = FALSE] %>%
    # rowSums() %>%
    apply(1, max) %>%
    {
      . * (-1)
    } %>%
    order()
  if (length(orders) < top_AEs) {
    AE_names <- rownames(x$contin_table)[filter_indi]
  } else {
    AE_names <- rownames(x$contin_table)[filter_indi][orders][1:top_AEs]
  }

  order_num_signal_per_drug <- (post_prob_matrix > 0.95) %>%
    colSums() %>%
    order()
  drug_names <- colnames(x$contin_table)[order_num_signal_per_drug] %>%
    rev()

  ordered_drug_names <- colnames(x$contin_table)[order_num_signal_per_drug]
  drug_names <- drug_names[1:top_drugs]
  if (!is.null(AEs)) {
    AE_names <- AEs
  }
  if (!is.null(drugs)) {
    drug_names <- drugs
  }
  dat_plot <- x$posterior_draws %>%
    # posterior::draws_of() %>%
    data.table::as.data.table() %>%
    data.table::setnames(
      old = c("V2", "V3", "V1", "value"),
      new = c("AE", "drug", "draw_idx", "post_draws")
    )

  dat_plot <- subset(dat_plot, (dat_plot$AE %in% AE_names) &
    (dat_plot$drug %in% drug_names)) %>%
    .left_join_base(counts_long, by = c("AE", "drug"))
  data.table::setDT(dat_plot)
  group_vars <- c("AE", "drug")
  measure_vars <- c("N", "E", "post_draws")
  dat_plot <- dat_plot[dat_plot$N > N_threshold, ]
  if (log_scale == TRUE) {
    dat_plot$post_draws[dat_plot$post_draws == 0] <- 1e-10
    dat_plot$post_draws <- log(dat_plot$post_draws)
    q05_cutoff <- log(1.01)
    xlab_text <- "Log signal strength (posterior median and 90% equi-tailed credible intervals)"
    vline_x <- 0
  } else {
    q05_cutoff <- 1.01
    xlab_text <- "Signal strength (posterior median and 90% equi-tailed credible intervals)"
    vline_x <- 1
  }


  dat_plot$AE <- (dat_plot$AE %>% .capitalize_words()) %>%
    factor(levels = AE_names %>% .capitalize_words() %>% rev())
  dat_plot$drug <- (dat_plot$drug %>% .capitalize_words()) %>%
    factor(levels = ordered_drug_names %>% .capitalize_words() %>% rev() %>%
      {
        c(setdiff(., "Other_drugs"), "Other_drugs")
      })


  q05_table <- dat_plot[, list(q05 = stats::quantile(.SD$post_draws, 0.05)), by = group_vars]

  q05_table <- q05_table[q05_table$q05 > q05_cutoff, ]
  dat_plot <- dat_plot[q05_table, on = group_vars]


  pl_summary <- dat_plot[, list(
    N = data.table::first(.SD$N),
    E = data.table::first(.SD$E),
    max_post_draws = max(.SD$post_draws),
    q95 = stats::quantile(.SD$post_draws, probs = c(0.95))
  ), by = group_vars, .SDcols = measure_vars]

  # Adding new columns using :=
  pl_summary[, `:=`(
    count_label = paste0("   N=", .SD$N),
    E_label = paste0("   E=", .SD$E)
  ), .SDcols = c("N", "E")]




  x_limit <- max(pl_summary$q95)

  pl <- dat_plot %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data[["post_draws"]],
        y = .data[["AE"]],
        group = .data[["drug"]],
        color = .data[["drug"]]
      )
    ) +
    ggdist::stat_pointinterval(
      position = ggplot2::position_dodge(0.9),
      .width = c(0.9),
      point_interval = "median_qi"
    ) +
    ggplot2::scale_x_continuous(
      # trans = "log1p",
      limits = c(NA, x_limit * x_lim_scalar)
    ) + # , limits = c(0, 13.5), breaks = c(0,1,3,5,7,10)) +
    ggplot2::geom_vline(
      xintercept = vline_x,
      color = "red",
      linetype = 2
    ) +
    ggplot2::geom_text(
      data = pl_summary,
      ggplot2::aes(
        x = max(.data[["q95"]]), y = .data[["AE"]],
        label = .data[["count_label"]], group = .data[["drug"]],
        color = .data[["drug"]]
      ), # Position the text on the right side
      position = ggplot2::position_dodge(0.9),
      hjust = 0,
      size = text_size,
      show.legend = FALSE,
      inherit.aes = FALSE # Prevents the legend issue
    ) +
    ggplot2::geom_text(
      data = pl_summary,
      ggplot2::aes(
        x = max(.data[["q95"]]) + text_shift, y = .data[["AE"]],
        label = .data[["E_label"]], group = .data[["drug"]],
        color = .data[["drug"]]
      ),
      position = ggplot2::position_dodge(0.9),
      hjust = 0,
      size = text_size,
      show.legend = FALSE,
      inherit.aes = FALSE # Prevents the legend issue
    ) +
    ggplot2::xlab(xlab_text) +
    wacolors::scale_color_wa_d(palette = "rainier") +
    ggplot2::theme_bw() +
    ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE))
  pl
}


.capitalize_words <- function(x) {
  vapply(x, function(s) {
    words <- strsplit(tolower(s), "\\s+")[[1]]
    capitalized <- paste(toupper(substring(words, 1, 1)),
      substring(words, 2),
      sep = ""
    )
    paste(capitalized, collapse = " ")
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
}


#' Generate a heatmap plot visualizing posterior probabilities for selected
#' drugs and adverse events
#'
#' @description
#' This function generates a heatmap to visualize the posterior probabilities
#' of being a signal for selected AEs and drugs.
#'
#'
#' @param x a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param num_top_AEs number of most significant AEs appearing in the plot. Default to
#' 10.
#' @param num_top_drugs number of most significant drugs appearing in the plot. Default
#' to 7.
#' @param specified_AEs a vector of AE names that are specified to appear in the plot. If
#' a vector of AEs is given, argument num_top_AEs will be ignored.
#' @param specified_drugs a vector of drug names that are specified to appear in the plot.
#' If a vector of drugs is given, argument num_top_drugs will be ignored.
#' @param cutoff_signal numeric. Threshold for signal detection. An AE-drug
#' combination is classified as a detected signal if its 5th posterior
#' percentile exceeds this threshold.
#'
#' @return
#' a ggplot2 object.
#'
#' @export
#'
#'
#' @examples
#'
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = 1000
#' )
#'
#'
#' heatmap_pvEBayes(
#'   x = fit,
#'   cutoff = 1.001
#' )
#'
heatmap_pvEBayes <- function(x,
                             num_top_AEs = 10,
                             num_top_drugs = 8,
                             specified_AEs = NULL,
                             specified_drugs = NULL,
                             cutoff_signal = NULL) {
  top_drugs <- num_top_drugs
  top_AEs <- num_top_AEs
  AEs <- specified_AEs
  drugs <- specified_drugs
  cutoff <- cutoff_signal
  if (top_drugs > (ncol(x$E) - 1)) {
    top_drugs <- (ncol(x$E) - 1)
  }
  if (is.null(cutoff)) {
    cutoff <- 1.001
  }



  stopifnot(is.pvEBayes(x))
  if (is.null(x$posterior_draws)) {
    x <- x %>% posterior_draws
  }
  counts_long <- x$contin_table %>%
    as.data.frame() %>%
    .rownames_to_column(var = "AE") %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars = "AE", variable.name = "drug", value.name = "N")
  Es_long <- x$E


  Es_long <- Es_long %>%
    round(2) %>%
    as.data.frame() %>%
    .rownames_to_column(var = "AE") %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars = "AE", variable.name = "drug", value.name = "E")

  counts_long <- counts_long %>%
    .left_join_base(Es_long, by = c("AE", "drug"))


  post_prob_matrix <- x$posterior_draws %>%
    {
      . > cutoff
    } %>%
    apply(c(2, 3), mean)

  RMSE1 <- x$posterior_draws %>%
    # posterior::draws_of() %>%
    {
      (. - 1)^2
    } %>%
    apply(c(2, 3), mean)

  filter_indi <- .check_AEs(post_prob_matrix, x$contin_table, 4)
  RMSE1[x$contin_table <= 4] <- 0
  orders <- RMSE1[filter_indi, , drop = FALSE] %>%
    # rowSums() %>%
    apply(1, max) %>%
    {
      . * (-1)
    } %>%
    order()
  if (length(orders) < top_AEs) {
    AE_names <- rownames(x$contin_table)[filter_indi]
  } else {
    AE_names <- rownames(x$contin_table)[filter_indi][orders][1:top_AEs]
  }


  order_num_signal_per_drug <- (post_prob_matrix > 0.95) %>%
    colSums() %>%
    order()
  drug_names <- colnames(x$contin_table)[order_num_signal_per_drug] %>%
    rev()

  ordered_drug_names <- colnames(x$contin_table)[order_num_signal_per_drug]
  drug_names <- drug_names[1:top_drugs]
  if (!is.null(AEs)) {
    AE_names <- AEs
  }
  if (!is.null(drugs)) {
    drug_names <- drugs
  }

  dat_plot <- post_prob_matrix %>%
    as.data.frame() %>%
    .rownames_to_column(var = "AE") %>%
    data.table::as.data.table() %>%
    data.table::melt(
      id.vars = "AE",
      variable.name = "drug",
      value.name = "post_prob"
    )


  dat_plot <- subset(
    dat_plot,
    (toupper(dat_plot$AE) %in% toupper(AE_names)) &
      (toupper(dat_plot$drug) %in% toupper(drug_names))
  ) %>%
    .left_join_base(counts_long, by = c("AE", "drug"))

  dat_plot$AE <- (dat_plot$AE %>% .capitalize_words()) %>%
    factor(levels = AE_names %>% .capitalize_words() %>% rev())
  dat_plot$drug <- (dat_plot$drug %>% .capitalize_words()) %>%
    factor(levels = ordered_drug_names %>% .capitalize_words() %>% rev())

  data.table::setDT(dat_plot)
  # Adding new columns using :=
  dat_plot[, `:=`(
    count_label = paste0(
      "N=", .SD$N, "; E=", .SD$E,
      "\n", "Post prob=", .SD$post_prob
    )
  ), .SDcols = c("N", "E", "post_prob")]



  pl <- dat_plot %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data[["drug"]],
        y = .data[["AE"]],
        fill = .data[["post_prob"]],
        label = .data[["count_label"]]
      )
    ) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradientn(colors = c("white", "blue")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1
      ),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "", y = "") +
    ggfittext::geom_fit_text(
      reflow = TRUE,
      contrast = TRUE,
      grow = TRUE
    )

  pl
}



#' Print method for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{print} method for objects of class
#' \code{pvEBayes}. It displays a concise description of the fitted model.
#'
#'
#' @param x a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#'
#' @param ... other input parameters. Currently unused.
#'
#'
#' @return
#' Invisibly returns the input `pvEBayes` object.
#'
#' @export
#'
#' @examples
#'
#' obj <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.5, n_posterior_draws = 10000
#' )
#'
#' print(obj)
print.pvEBayes <- function(x, ...) {
  stopifnot(is.pvEBayes(x))
  top_text <- "Object of class 'pvEBayes'"
  if (x$model == "general-gamma") {
    model_info_text <- glue::glue(
      "General-gamma model with hyperparameter alpha = ",
      x$alpha, ".\n", "Estimated prior is a mixture of ",
      length(x$r),
      " gamma distributions."
    )
  } else if (x$model == "K-gamma") {
    model_info_text <- glue::glue(
      "K-gamma model with number of gamma mixture components K = ",
      length(x$r), "."
    )
  } else if (x$model == "GPS") {
    model_info_text <- glue::glue("GPS (2-gamma) model is fitted")
  } else if (x$model == "KM") {
    model_info_text <- glue::glue("KM model is fitted")
  } else {
    model_info_text <- glue::glue(
      "efron model is fitted with hyperparameters",
      " (p = ", x$p, ", c0 = ", x$c0, ")."
    )
  }

  run_time_txt1 <- glue::glue(
    "Running time of the ", x$model,
    " model fitting: ",
    x$fit_time %>%
      as.numeric(units = "secs") %>%
      round(4),
    " seconds."
  )
  if (!is.null(x$draws_time)) {
    run_time_txt2 <- glue::glue(
      "Running time for posterior draws \n (",
      x$n_posterior_draws,
      " signal strength posterior draws per AE-drug pair):",
      x$draws_time %>%
        as.numeric(units = "secs") %>%
        round(4),
      " seconds."
    )
  } else {
    run_time_txt2 <- "No posterior draws were generated."
  }
  run_time_txt <- glue::glue(
    "{run_time_txt1}
    {run_time_txt2}
    "
  )
  msg <- glue::glue(
    "{top_text}

    {model_info_text}

    {run_time_txt}

    Extract estimated prior parameters, discovered signals
    and signal strength posterior draws using `summary()`.
    "
  )
  message(msg)
  invisible(x)
}



#' Summary method for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{summary} method for objects of class
#' \code{pvEBayes}. It provides a detailed summary of the fitted model.
#'
#'
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#'
#' @param return A character string specifying which component the summary
#'  function should return.Valid options include: "prior parameters",
#' "likelihood", "detected signal" and "posterior draws". If set to NULL
#' (default), all components will be returned in a list.
#'
#' @param ... other input parameters. Currently unused.
#'
#' @return
#' a list including estimated prior parameters, log_marginal_likelihood,
#' indicator matrix of detected signal and posterior_draws for each AE-drug
#' pair.
#'
#' @export
#'
#' @examples
#'
#' obj <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.5, n_posterior_draws = 10000
#' )
#'
#' summary(obj)
summary.pvEBayes <- function(object, return = NULL,...) {
  stopifnot(is.pvEBayes(object))
  if (is.null(object$posterior_draws)) {
    object <- posterior_draws(object, n_posterior_draws = 1000)
  }
  model <- object$model
  if (model == "KM") {
    estimated_prior <- list(
      g = object$g,
      grid = object$grid
    )
  } else if (model == "efron") {
    estimated_prior <- list(
      a = object$a,
      g = object$g,
      grid = object$grid
    )
  } else {
    estimated_prior <- list(
      omega = object$omega,
      r = object$r,
      h = object$h
    )
  }
  if(!is.null(return)){
    if(return == "prior parameters"){
      estimated_prior
    }else if(return == "likelihood"){
      object$loglik
    }else if(return == "detected signal"){
      (.get_posterior_prob(object) >= 0.95)
    }else if(return ==  "posterior draws"){
      object$posterior_draws
    }else{
      stop("Please provide a valid return argument.")
    }
  }else{
    res <- list(
      estimated_prior = estimated_prior,
      log_marginal_likelihood = object$loglik,
      detected_signal = (.get_posterior_prob(object) >= 0.95),
      posterior_draws = object$posterior_draws
    )
    res
  }


}



#' Plotting method for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{plot} method for objects of class
#' \code{pvEBayes}.
#'
#'
#' @param x a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param type character string determining the type of plot to show.
#' Available choices are `"eyeplot"` which calls \link{eyeplot_pvEBayes} and
#' `"heatmap"` which calls \link{heatmap_pvEBayes}, with the additional arguments
#' supplied in ...
#' @param ... additional arguments passed to heatmap_pvEBayes or eyeplot_pvEBayes.
#'
#' @return
#'
#' A \link[ggplot2]{ggplot} object.
#' @export
#'
#' @examples
#'
#' obj <- pvEBayes(statin2025_44, model = "general-gamma", alpha = 0.5)
#' plot(obj, type = "eyeplot")
plot.pvEBayes <- function(x, type = "eyeplot", ...) {
  if (!is.pvEBayes(x)) {
    stop("x must be a 'pvEBayes' object.")
  }

  stopifnot(
    type %in% c("heatmap", "eyeplot")
  )

  out <- if (type == "heatmap") {
    heatmap_pvEBayes(x, ...)
  } else if (type == "eyeplot") {
    eyeplot_pvEBayes(x, ...)
  }

  out
}


#' Extract log marginal likelihood for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{logLik} method for objects of class
#' \code{pvEBayes}. It extracts the log marginal likelihood from a fitted
#' model.
#'
#'
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param ... other input parameters. Currently unused.
#'
#' @return
#' returns log marginal likelihood of a pvEBayes object.
#'
#' @export
#'
#' @examples
#'
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = NULL
#' )
#'
#' logLik(fit)
#'
logLik.pvEBayes <- function(object, ...) {
  stopifnot(is.pvEBayes(object))

  res <- object$loglik
  return(res)
}




#' Obtain Akaike Information Criterion (AIC) for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{AIC} method for objects of class
#' \code{pvEBayes}. It extracts the Akaike Information Criterion (AIC)
#' from a fitted model.
#'
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param ... other input parameters. Currently unused.
#' @param k numeric, the penalty per parameter to be used; the default k = 2
#' is the classical AIC.
#'
#' @return
#' numeric, AIC score for the resulting model.
#'
#' @export
#'
#' @examples
#'
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = NULL
#' )
#'
#' AIC_score <- AIC(fit)
#'
AIC.pvEBayes <- function(object, ..., k = 2) {
  stopifnot(is.pvEBayes(object))
  model <- object$model
  if (model == "KM") {
    penalty <- length(object$g)
  } else if (model == "efron") {
    penalty <- object$df
  } else {
    penalty <- length(object$r) * 3
  }


  AIC_score <- penalty * 2 - object$loglik * k
  return(AIC_score)
}


#' Obtain Bayesian Information Criterion (BIC) for a pvEBayes object
#'
#' @description
#' This function defines the S3 \code{BIC} method for objects of class
#' \code{pvEBayes}. It extracts the Bayesian Information Criterion (BIC)
#' from a fitted model.
#'
#'
#' @param object a \code{pvEBayes} object, which is the output of the function
#' \link{pvEBayes} or \link{pvEBayes_tune}.
#' @param ... other input parameters. Currently unused.
#'
#' @return
#'
#' numeric, BIC score for the resulting model.
#'
#' @importFrom stats BIC
#' @export
#'
#' @examples
#'
#' fit <- pvEBayes(
#'   contin_table = statin2025_44, model = "general-gamma",
#'   alpha = 0.3, n_posterior_draws = NULL
#' )
#'
#' BIC_score <- BIC(fit)
#'
BIC.pvEBayes <- function(object, ...) {
  stopifnot(is.pvEBayes(object))
  model <- object$model
  n <- object$contin_table %>%
    dim() %>%
    prod()

  if (model == "KM") {
    penalty <- length(object$g)
  } else if (model == "efron") {
    penalty <- object$df
  } else {
    penalty <- length(object$r) * 3
  }


  BIC_score <- penalty * log(n) - object$loglik * 2
  return(BIC_score)
}
