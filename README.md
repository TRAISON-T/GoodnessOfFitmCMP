[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) <!-- badges: start --> [![R-CMD-check](https://github.com/TRAISON-T/GoodnessOfFitmCMP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TRAISON-T/GoodnessOfFitmCMP/actions/workflows/R-CMD-check.yaml) <!-- badges: end --> <!-- badges: start --> [![Codecov test coverage](https://codecov.io/gh/TRAISON-T/GoodnessOfFitmCMP/graph/badge.svg)](https://app.codecov.io/gh/TRAISON-T/GoodnessOfFitmCMP) <!-- badges: end --> [![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# GoodnessOfFitmCMP

GoodnessOfFitmCMP is an R package that provides goodness-of-fit tests for the **Conway-Maxwell Poisson (COM-Poisson)** distribution, specifically its mean-reparametrized form (mCMP). The package implements novel test statistics based on **Stein's characterization**, which offers a computationally efficient way to handle the normalizing constant inherent in the COM-Poisson probability mass function.

## Overview

The COM-Poisson distribution is a flexible model for count data that can accommodate **equi-dispersion** (like Poisson), **under-dispersion**, and **over-dispersion**. Traditional goodness-of-fit tests for this distribution can be computationally challenging. This package addresses that by using Stein's characterization to define modified test statistics:

-   **Modified Cramér-von Mises** (CV_M)
-   **Modified Anderson-Darling** (AD_M)
-   **Probability Distance** (PD)

The package also includes implementations of the standard Cramér-von Mises (CV), Anderson-Darling (AD), and Chi-square (χ²) tests for comparison. P-values for these tests are computed using a **bootstrap methodology**.

## Installation

You can install the `GoodnessOfFitmCMP` package directly from GitHub (if you host it there) or locally using `devtools`.

First, ensure you have the `devtools` package installed:

``` r
install.packages("devtools")
```

### Install from GitHub

``` r
devtools::install_github("TRAISON-T/GoodnessOfFitmCMP")
```

## Key Functions

-   `dcmp(x, mu, phi, max_terms, tol)`: Computes the **probability mass function (PMF)** of the mCMP distribution.
-   `simulate_mcmp_sample(n, mu, phi, max_terms, tol)`: Generates **random samples** from the mCMP distribution using inverse transform sampling.
-   `goodness_of_fit_mCMP(x, B, verbose)`: The **main function** to perform the goodness-of-fit tests.
-   `estimate_mu_phi(x)`: Estimates the mCMP parameters μ and φ.
-   `S_objective(par, x)`: The objective function minimized for parameter estimation.
-   `estimate_rho_mCMP(x, mu, phi, k)`: Estimates the **Stein-based probability** ρ_X(k).
-   `compute_hat_rho_vec(x, mu, phi, r)`: Computes vectorized Stein-based probabilities.
-   `compute_mod_statistics(sample, mu, phi)`: Calculates the various goodness-of-fit test statistics.

For detailed documentation on each function, use `?function_name` in R (e.g., `?dcmp`) after installing the package.

## Usage

``` r
# Load the package
library(GoodnessOfFitmCMP)

# Example data 
x <- c(12, 8, 9, 11, 9, 10, 12, 9, 12, 4, 11, 8, 12, 12, 9, 8, 9, 10, 15, 11, 14, 11, 9, 7, 15, 11, 13, 9, 15, 8, 7, 11)

# Perform the goodness-of-fit test
# Use a smaller B for quick examples; B=1000 or more is typically recommended for real analysis.
gof_result <- goodness_of_fit_mCMP(x, B = 100, verbose = FALSE) # verbose=FALSE to suppress bootstrap messages

# Print the results
print(gof_result)

# Access estimated parameters
cat("Estimated Mu:", gof_result$estimated_mu, "\n")
cat("Estimated Phi:", gof_result$estimated_phi, "\n")

# Access p-values
print(gof_result$p_value)

# Access observed test statistics
print(gof_result$obs_stats)
```

## Theoretical Background

This package is based on the research presented in the article "Goodness-of-Fit Tests for COM-Poisson Distribution Using Stein's Characterization" by Traison T. and V.S. Vaidyanathan. The core idea is to leverage **Stein's characterization** to represent the COM-Poisson PMF in a form that avoids the computationally intensive normalizing constant. This allows for the development of new, more efficient goodness-of-fit test statistics. The parameters are estimated by minimizing a probability distance criterion, and the significance of the tests is assessed via a **bootstrap resampling procedure**.

## License

This package is released under the MIT License. See the `LICENSE` file for more details.
