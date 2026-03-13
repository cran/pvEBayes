---
output:
  html_document: default
  pdf_document: default
---

# R package `pvEBayes`

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pvEBayes)](https://CRAN.R-project.org/package=pvEBayes)
[![R-CMD-check](https://github.com/YihaoTancn/pvEBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YihaoTancn/pvEBayes/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/YihaoTancn/pvEBayes/graph/badge.svg)](https://app.codecov.io/gh/YihaoTancn/pvEBayes)
[![CodeFactor](https://www.codefactor.io/repository/github/yihaotancn/pvebayes/badge)](https://www.codefactor.io/repository/github/yihaotancn/pvebayes)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!--[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/746_status.svg)](https://github.com/ropensci/software-review/issues/746)-->
<!-- badges: end -->

`pvEBayes` is an R package that implements a suite of nonparametric empirical
Bayes methods for pharmacovigilance, including the Gamma-Poisson Shrinker (GPS),
K-gamma, general-gamma, Koenker-Mizera (KM), and Efron models. It provides tools
for fitting these models to the spontaneous reporting system (SRS) frequency 
tables, extracting summaries, performing hyperparameter tuning, and generating 
graphical summaries (eye plots and heatmaps) for signal detection and signal 
strength estimation.

**Spontaneous Reporting System (SRS) Table**: An drug safety SRS dataset 
catalogs AE reports on *I* AE rows across *J* drug columns. Let ${N_{ij}}$ 
denote the number of reported cases for the *i*-th AE and the *j*-th drug, 
where ${i = 1,..., I}$ and ${j = 1,..., J}$. 

**Empirical Bayes modeling for disproportionality analysis**: 

  - Model each AE-drug count as 
  $N_{ij} \sim \text{Poisson}(\lambda_{ij} E_{ij})$, $N_{ij} = 0, 1, 2, \dots$
  
  - $E_{ij}$: expected baseline value assuming no AE-drug association.
  
  - $\lambda_{ij} > 0$: relative reporting ratio / signal strength for the 
  (i,j)-th AE-drug pair (multiplicative deviation from the null baseline value).


**From signal detection to signal strength $\lambda$ estimation**

  - Traditional disproportionality analysis emphasizes **signal detection**: 
  identify AE-drug pairs with observed counts substantially larger than 
  its null value, i.e., $\lambda_{ij} > 1$.
  
  - Tan et al. (*Stat. in Med.*, 2025) extend this to 
  **signal strength estimation**: estimate $\{\lambda_{ij}\}$ and quantify 
  uncertainty via flexible nonparametric empirical Bayes posterior distribution.
  
  - Signal estimation helps distinguish AE-drug pairs that look identical under 
  a binary signal detection framework (e.g., $\lambda=1.5$ vs $\lambda=4.0$), 
  which can have different clinical implications.

**Methods implemented in `pvEBayes` (differ by prior assumptions):**

  - Gamma-Poisson Shrinker (GPS)
  
  - Koenker–Mizera method (KM)
  
  - Efron’s nonparametric empirical Bayes approach
  
  - K-gamma method
  
  - General-gamma method
  
**Why nonparametric priors?**

  - GPS uses a 2 gamma mixture prior motivated by a signal/non-signal structure.
  
  - Real-world signal strengths $\lambda_{ij}$ can be heterogeneous, and the 
  underlying (prior) distribution over $\lambda_{ij}$ may be multi-modal with
  multiple distinct peaks, making simple parametric priors hard to justify.
  
  - Nonparametric empirical Bayes methods (KM, Efron, general-gamma) address 
  this challenge by utilizing a flexible prior with a general mixture form and 
  estimating the prior distribution in a data-driven way. 

**Implementation highlights:**

  - Provides a fully open-source KM implementation using `CVXR` (avoids reliance 
  on the commercial Mosek solver used by `REBayes`).
  
  - Adapts Efron’s approach from `deconvolveR` to support the exposure/offset 
  $E_{ij}$ in the Poisson model (not supported in the original implementation).
  
  - Implements the bi-level Expectation Conditional Maximization (ECM) algorithm 
  from Tan et al. (*Stat. in Med.*, 2025) for prior estimation for 
  gamma-mixture prior based models (GPS, K-gamma, general-gamma).
  
For a detailed methodological description, see Tan et al. 
(*Stat. in Med.*, 2025).


## Installation

The stable version of `pvEBayes` can be installed from CRAN:

```
install.packages("pvEBayes")
```

The development version is available from GitHub:

```
# if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("YihaoTancn/pvEBayes")
```

## Quick Example

Here is a minimal example analyzing the built-in FDA statin44 dataset with
general-gamma model:

```

library(pvEBayes)

# Load the statin44 contingency table of 44 AEs for 6 statins
data("statin2025_44")

# Fit a general-gamma model with a specified alpha
fit <- pvEBayes(
  contin_table      = statin2025_44,
  model             = "general-gamma",
  alpha             = 0.3,
  n_posterior_draws = 1000
)

# Print out a concise description of the fitted model
fit

# Obtain a logical matrix for the detected signal
summary(fit, return = "detected signal")

# Visualize posterior distributions for top AE-drug pairs
plot(fit, type = "eyeplot")

```

For a more detailed illustration, please see 'Vignette'.

## License

`pvEBayes` is released under the GPL-3 license. See 'LICENSE.md' for details.


## Code of Conduct
  
Please note that the `pvEBayes` project is released with a 
[Contributor Code of Conduct](https://yihaotancn.github.io/pvEBayes/CODE_OF_CONDUCT.html). 
By contributing to this project, you agree to abide by its terms.

## References

Tan Y, Markatou M and Chakraborty S. Flexible Empirical Bayesian Approaches
to Pharmacovigilance for Simultaneous Signal Detection and Signal Strength
Estimation in Spontaneous Reporting Systems Data.
*Statistics in Medicine*. 2025; 44: 18-19,
https://doi.org/10.1002/sim.70195.

Tan Y, Markatou M and Chakraborty S. pvEBayes: An R Package for Empirical Bayes 
Methods in Pharmacovigilance. *arXiv*:2512.01057 (stat.AP). 
https://doi.org/10.48550/arXiv.2512.01057

Koenker R, Mizera I. Convex Optimization, Shape Constraints, Compound
Decisions, and Empirical Bayes Rules. *Journal of the American
Statistical Association* 2014; 109(506): 674–685,
https://doi.org/10.1080/01621459.2013.869224

Efron B. Empirical Bayes Deconvolution Estimates. *Biometrika* 2016;
103(1); 1-20, https://doi.org/10.1093/biomet/asv068

DuMouchel W. Bayesian data mining in large frequency tables, with an
application to the FDA spontaneous reporting system.
*The American Statistician*. 1999; 1;53(3):177-90.
