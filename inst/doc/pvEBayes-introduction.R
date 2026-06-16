## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(pvEBayes)
library(ggplot2)

set.seed(1)
# load the SRS data
data("statin2025_44")

# show the first 6 rows
head(statin2025_44)

## ----include=TRUE-------------------------------------------------------------
gg_given_alpha <- pvEBayes(statin2025_44,
  model = "general-gamma",
  alpha = 0.5,
  maxi = 200,
  tol_ecm = 1e-4  #default value
)

gg_given_alpha$iter

## ----include=TRUE-------------------------------------------------------------


gg_given_alpha2 <- pvEBayes(statin2025_44,
  model = "general-gamma",
  alpha = 0.5,
  maxi = 200,
  tol_ecm = 1e-8  #smaller tolerance for convergence 
)

gg_given_alpha2$iter

## ----fit-gg, include=TRUE-----------------------------------------------------
summary(gg_given_alpha)

gg_given_alpha_detected_signal <- summary(gg_given_alpha,
  return = "detected signal"
)
sum(gg_given_alpha_detected_signal)

## ----get-post, include=TRUE---------------------------------------------------
gg_customize_detected_signal <- get_posterior_prob(gg_given_alpha,
                                                     cutoff_signal = 1.01) > 0.99
sum(gg_customize_detected_signal)

## ----tune, include=TRUE-------------------------------------------------------
AIC(gg_given_alpha)

BIC(gg_given_alpha)

## ----auto-tune, include=TRUE, eval=TRUE---------------------------------------
gg_tune_statin44 <- pvEBayes_tune(statin2025_44,
  model = "general-gamma",
  alpha_vec = c(0, 0.1, 0.3, 0.5, 0.7, 0.9),
  use_AIC = TRUE
)

## ----e-auto-tune, include=TRUE, eval=TRUE-------------------------------------

e_tune_statin44 <- pvEBayes_tune(statin2025_44,
  model = "efron",
  p_vec = c(40, 60, 80),
  c0 = c(0.001, 0.01, 0.1),
  use_AIC = TRUE
)

e_tune_statin44

## ----heatmap, include = TRUE, fig.height=6, fig.width=8-----------------------
heatmap_gg_tune_statin44 <- plot(gg_tune_statin44,
  type = "heatmap",
  num_top_AEs = 10,
  cutoff_signal = 1.001
)

heatmap_gg_tune_statin44 +
  theme(
    legend.position = "top"
  )

## ----eyeplot, include = TRUE, fig.height=6, fig.width=8-----------------------
eyeplot_gg_tune_statin44 <- plot(gg_tune_statin44,
  type = "eyeplot",
  num_top_AEs = 8,
  N_threshold = 1,
  log_scale = FALSE,
  text_shift = 2.3,
  text_size = 3,
  x_lim_scalar = 1.2
)

eyeplot_gg_tune_statin44 +
  theme(
    legend.position = "top"
  )

