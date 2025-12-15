
# R package `pvEBayes`

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pvEBayes)](https://CRAN.R-project.org/package=pvEBayes)
[![R-CMD-check](https://github.com/YihaoTancn/pvEBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YihaoTancn/pvEBayes/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/YihaoTancn/pvEBayes/graph/badge.svg)](https://app.codecov.io/gh/YihaoTancn/pvEBayes)
[![CodeFactor](https://www.codefactor.io/repository/github/yihaotancn/pvebayes/badge)](https://www.codefactor.io/repository/github/yihaotancn/pvebayes)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`pvEBayes` is an R package that implements a suite of nonparametric empirical
Bayes methods for pharmacovigilance, including Gamma-Poisson Shrinker (GPS),
K-gamma, general-gamma, Koenker-Mizera (KM), and Efron models. It provides tools
for fitting these models to the spontaneous reporting system (SRS) frequency 
tables, extracting summaries, performing hyperparameter tuning, and generating 
graphical summaries (eye plots and heatmaps) for signal detection and 
estimation.

**Spontaneous Reporting System (SRS) Table**: An drug safety SRS dataset 
catalogs AE reports on *I* AE rows across *J* drug columns. Let ${N_{ij}}$ 
denote the number of reported cases for the *i*-th AE and the *j*-th drug, 
where ${i = 1,..., I}$ and ${j = 1,..., J}$. 

**Empirical Bayes modeling for disproportionality analysis**: We assume that for 
each AE-drug pair, $N_{ij} \sim \text{Poisson}(\lambda_{ij} E_{ij})$, where 
${E_{ij}}$ is the expected baseline value, measuring the expected count of the 
AE-drug pair when there is no association between*i*-th AE and *j*-th drug. 
The parameter ${\lambda_{ij} \geq 0}$ represents the relative reporting ratio, 
the signal strength, for the ${(i, j)}$-th pair, measuring the ratio of the 
actual expected count arising due to dependence on the null baseline expected 
count. Current disproportionality analysis mainly focuses on *signal detection* 
which seeks to determine whether the observation $N_{ij}$ is substantially 
greater than the corresponding null baseline $E_{ij}$. Under the Poisson model, 
that is to say, its signal strength $\lambda_{ij}$ is significantly greater 
than 1.

In addition to *signal detection*, Tan et al. (*Stat. in Med.*, 2025) broaden 
the role of disproportionality to *signal estimation*. The use of the flexible 
non-parametric empirical Bayes models enables more nuanced empirical Bayes 
posterior inference (parameter estimation and uncertainty quantification) on 
signal strength parameter $\{ \lambda_{ij} \}$. This allows researchers to 
distinguish AE-drug pairs that would appear similar under a binary signal 
detection framework. For example, the AE-drug pairs with signal strengths of 
1.5 and 4.0 could both be significantly greater than 1 and detected as a 
signal. Such differences in signal strength may have distinct implications in 
medical and clinical contexts.

The methods included in `pvEBayes` differ by their assumptions on the prior 
distribution. Implemented methods include the Gamma-Poisson Shrinker (GPS),
Koenker-Mizera (KM) method, Efron’s nonparametric empirical Bayes approach,
the K-gamma model, and the general-gamma model. The selection of the prior 
distribution is critical in Bayesian analysis. The GPS model uses a gamma 
mixture prior by assuming the signal/non-signal structure in SRS data. 
However, in real-world setting, signal strengths ${(\lambda_{ij})}$ are 
often heterogeneous and thus follows a multi-modal distribution, making it
difficult to assume a parametric prior. Non-parametric empirical Bayes models 
(KM, Efron, K-gamma and general-gamma) address this challenge by utilizing a 
flexible prior with a general mixture form and estimating the prior 
distribution in a data-driven way.    


**Implementations**: The KM method has an existing implementation in the 
`REBayes` package, but it relies on Mosek, a commercial convex optimization 
solver, which may limit accessibility due to licensing issues. The `pvEBayes` 
package provides an alternative fully open-source implementation of the KM 
method using `CVXR`. Efron’s method also has a general nonparametric empirical 
Bayes implementation in the `deconvolveR` package; however, that
implementation does not support an exposure or offset parameter in the
Poisson model, which corresponds to the expected null value ${E_{ij}}$. In
`pvEBayes`, the implementation of Efron's method is adapted and
modified from `deconvolveR` to support ${E_{ij}}$ in the Poisson model.

In addition, this package implements the novel bi-level Expectation
Conditional Maximization (ECM) algorithm proposed by Tan et al. (2025) for
efficient parameter estimation in gamma mixture prior based models mentioned
above.



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
