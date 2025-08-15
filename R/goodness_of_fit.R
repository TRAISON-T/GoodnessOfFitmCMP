# GoodnessOfFitmCMP/R/goodness_of_fit.R
#' @title Goodness-of-Fit Tests for the mCMP Distribution
#'
#' @description Performs goodness-of-fit tests for the mean-reparametrized Conway-Maxwell Poisson (mCMP) distribution using bootstrap methodology.
#'   The mCMP distribution is a two-parameter generalization of the Poisson distribution capable of modeling
#'   equi-, under-, and over-dispersed count data. Assessing the adequacy of fit is crucial in data analysis.
#'   This function implements proposed modified goodness-of-fit tests based on Stein's characterization:
#'   the modified Cramér-von Mises (\eqn{CV_M}), modified Anderson-Darling (\eqn{AD_M}),
#'   and Probability Distance (\eqn{PD}) tests. These tests overcome the computational complexity
#'   associated with the normalizing constant in the mCMP probability mass function.
#'   For comparison, it also includes the traditional Cramér-von Mises (CV),
#'   Anderson-Darling (AD), and Chi-square (\eqn{\chi^2}) tests.
#'
#'   The p-values for these test statistics are computed using a bootstrap procedure. This is necessary
#'   because the exact distributions of these test statistics under the null hypothesis
#'   (that the data comes from an mCMP distribution) are generally difficult to obtain.
#'   The parameters \eqn{\mu} and \eqn{\phi} of the mCMP distribution are estimated internally
#'   by minimizing the Probability Distance (\eqn{PD}) objective function, \eqn{S(\mu,\phi)}.
#'
#' @param x A numeric vector of observed sample data (non-negative integers).
#' @param B An integer for the number of bootstrap runs. Default is 1000. A larger number provides
#'   more accurate p-values but takes longer to compute.
#' @param verbose A logical value indicating whether to print progress messages during the bootstrap iterations.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{estimated_mu}: The estimated location parameter \eqn{\mu} of the mCMP distribution.
#'     \item \code{estimated_phi}: The estimated dispersion parameter \eqn{\phi} of the mCMP distribution.
#'     \item \code{p_value}: A named numeric vector of p-values for each test statistic (CV, AD, chi_sq, CV_M, AD_M, PD).
#'       These p-values indicate the probability of observing a test statistic as extreme or more extreme
#'       than the observed one, assuming the null hypothesis is true.
#'     \item \code{obs_stats}: A named numeric vector of the observed test statistics calculated from the input data \code{x}.
#'   }
#' @examples
#' # Example data (from the paper's simulation study)
#' x <- c(12, 8, 9, 11, 9, 10, 12, 9, 12, 4, 11, 8, 12, 12, 9, 8, 9, 10,
#'        15, 11, 14, 11, 9, 7, 15, 11, 13, 9, 15, 8, 7, 11)
#'
#' # Run the goodness-of-fit test. Use a smaller B for quick examples;
#' # B=1000 or more is typically recommended for real analysis.
#' gof_result <- goodness_of_fit_mCMP(x, B = 100, verbose = FALSE)
#' print(gof_result)
goodness_of_fit_mCMP <- function(x, B = 1000, verbose = TRUE) {
  n <- length(x)
  r <- max(x)

  # Step 1: Estimate mu, phi using PD criterion
  est <- estimate_mu_phi(x)
  mu_hat <- est[1]
  phi_hat <- est[2]

  # Step 2: Compute observed stats
  obs_stats <- compute_mod_statistics(x, mu_hat, phi_hat)

  # Step 3: Bootstrap
  boot_stats <- matrix(NA, nrow = B, ncol = 6)
  colnames(boot_stats) <- c("CV", "AD","chi_sq","CV_M", "AD_M", "PD")

  for (b in 1:B) {
    x_boot <- simulate_mcmp_sample(n, mu_hat, phi_hat)
    est1 <- estimate_mu_phi(x_boot)
    mu_hat1 <- est1[1]

    phi_hat1 <- est1[2]
    boot_stats[b, ] <- compute_mod_statistics(x_boot, mu_hat1, phi_hat1)
    if (verbose && b %% (B / 10) == 0) cat("Bootstrap iteration:", b, "of", B, "\n")
  }

  # Step 4: Compute p-values
  pvals <- colMeans(sweep(boot_stats, 2, obs_stats, `>=`))

  return(list(
    estimated_mu = mu_hat,
    estimated_phi = phi_hat,
    p_value = pvals, obs_stats = obs_stats
  ))
}

#' @title Estimator of rho_X(k) (Internal Helper)
#'
#' @description This internal helper function computes an estimator of \eqn{\rho_X(k)}, which represents
#'   the probability mass function (PMF) of the mCMP distribution derived using Stein's characterization.
#'   Stein's characterization is a powerful tool that allows the PMF to be expressed in a form
#'   that avoids the computationally challenging normalizing constant present in the original mCMP PMF.
#'   The estimator \eqn{\hat{\rho}_X(k)} is defined as:
#'   \deqn{ \hat{\rho}_X(k) = \frac{1}{n}\sum_{i=1}^n \left( 1 - \left(\frac{\mu + \frac{e^\phi - 1}{2e^\phi}}{x_i + 1}\right)^{e^\phi} \right) I(x_i \ge k) }
#'   where \eqn{x_i} are the observed sample values, \eqn{n} is the sample size, and \eqn{I(x_i \ge k)}
#'   is the indicator function, which is 1 if \eqn{x_i \ge k} and 0 otherwise.
#'   By the Law of Large Numbers, \eqn{\hat{\rho}_X(k)} converges in probability to \eqn{\rho_X(k)}.
#'
#' @param x A numeric vector of sample data.
#' @param mu Numeric location parameter \eqn{\mu}.
#' @param phi Numeric dispersion parameter \eqn{\phi}.
#' @param k Integer for which \eqn{\rho_X(k)} is estimated.
#'
#' @return Numeric value of the estimator. Returns 0 if parameters lead to invalid calculations (e.g., \code{constant} is non-positive).
#' @examples
#' # Simulate some mCMP data
#' set.seed(123)
#' sample_data <- simulate_mcmp_sample(n = 50, mu = 5, phi = 0.1)
#'
#' # Estimate rho_X(k) for specific k values
#' # Estimate rho_X(0) (probability of observing 0)
#' rho_0 <- estimate_rho_mCMP(x = sample_data, mu = 5, phi = 0.1, k = 0)
#' cat("Estimated rho_X(0):", rho_0, "\n")
#'
#' # Estimate rho_X(1) (probability of observing 1)
#' rho_1 <- estimate_rho_mCMP(x = sample_data, mu = 5, phi = 0.1, k = 1)
#' cat("Estimated rho_X(1):", rho_1, "\n")
#'
#' # Note: For internal use, this function is typically called by compute_hat_rho_vec.
estimate_rho_mCMP <- function(x, mu, phi, k) {
  offset <- (exp(phi) - 1) / (2 * exp(phi))
  constant <- mu + offset
  if (constant <= 0) return(0)
  term <- 1 - (constant / (x + 1))^exp(phi)
  indicator <- as.numeric(x >= k)
  val <- mean(term * indicator)
  if(is.na(val)||val<0){val=0}
  return(val)
}

#' @title Vectorized rho_X for k = 0 to r (Internal Helper)
#'
#' @description This internal helper function efficiently computes the estimated Stein-based probabilities
#'   \eqn{\hat{\rho}_X(k)} for a range of \eqn{k} values, specifically from 0 up to the maximum observed
#'   value \eqn{r} in the sample. It uses the \code{\link{estimate_rho_mCMP}} function for each \eqn{k}.
#'   This vectorized approach is used to construct the estimated PMF based on Stein's characterization.
#'
#' @param x Numeric sample data.
#' @param mu Numeric location parameter \eqn{\mu}.
#' @param phi Numeric dispersion parameter \eqn{\phi}.
#' @param r Maximum value of \eqn{k} to compute \eqn{\hat{\rho}_X(k)} for. This is typically \code{max(x)}
#'   from the input sample.
#'
#' @return A numeric vector of \eqn{\hat{\rho}_X(k)} values for \eqn{k = 0, \ldots, r}.
#' @examples
#' # Simulate some mCMP data
#' set.seed(123)
#' sample_data <- simulate_mcmp_sample(n = 20, mu = 3, phi = -0.5)
#'
#' # Compute rho_X(k) for k from 0 to max(sample_data)
#' max_val <- max(sample_data)
#' hat_rho_values <- compute_hat_rho_vec(x = sample_data, mu = 3, phi = -0.5, r = max_val)
#' names(hat_rho_values) <- 0:max_val # Assign names for clarity
#' print(hat_rho_values)
#'
#' # These values represent the estimated PMF based on Stein's characterization.
compute_hat_rho_vec <- function(x, mu, phi, r) {
  sapply(0:r, function(k) estimate_rho_mCMP(x, mu, phi, k))
}

#' @title Compute Goodness-of-Fit Test Statistics (Internal Helper)
#'
#' @description This internal helper function calculates various goodness-of-fit test statistics
#'   to assess the fit of an mCMP distribution to observed sample data.
#'   It computes both traditional and modified statistics, as detailed in the article.
#'
#'   \strong{Traditional Statistics:}
#'   \itemize{
#'     \item \strong{Cramér-von Mises (CV):} Measures the squared difference between the empirical
#'       and theoretical cumulative distribution functions (CDFs), weighted by the theoretical PMF.
#'       \deqn{ CV = n\sum_{x = 0}^{r} {{{\left[ {{F_n}\left( x \right) - F\left( {x;\mu ,\phi } \right)} \right]}^2}p\left( {x;\mu ,\phi } \right)} }
#'     \item \strong{Anderson-Darling (AD):} Similar to CV, but gives more weight to the tails of the distribution.
#'       \deqn{ AD = n\sum_{x = 0}^{r} {\frac{{{{\left[ {{F_n}\left( x \right) - F\left( {x;\mu ,\phi } \right)} \right]}^2}p\left( {x;\mu ,\phi } \right)}}{{F\left( {x;\mu ,\phi } \right)\left[ {1 - F\left( {x;\mu ,\phi } \right)} \right]}}} }
#'     \item \strong{Chi-square (\eqn{\chi^2}):} Compares observed frequencies (\eqn{O_x}) with expected frequencies (\eqn{E_x})
#'       under the hypothesized distribution.
#'       \deqn{ \chi^2 = \sum_{x = 0}^{r} \frac{(O_x - E_x)^2}{E_x} }
#'   }
#'
#'   \strong{Modified Statistics (using Stein's characterization):}
#'   These statistics are proposed in the article to overcome computational challenges
#'   due to the normalizing constant in the original mCMP PMF. They use the Stein-based
#'   estimated PMF \eqn{\hat{\rho}_X(x)} and its corresponding CDF \eqn{F_{\hat{\rho}_X}(x)}.
#'   \itemize{
#'     \item \strong{Modified Cramér-von Mises (CV_M):}
#'       \deqn{ C{V_M} = n\sum_{x = 0}^{r} {{{\left[ {{F_n}\left( x \right) - {F_{{{\hat \rho }_x}}}\left( {x;\mu ,\phi } \right)} \right]}^2}{{\hat \rho }_{X}(x)}} }
#'     \item \strong{Modified Anderson-Darling (AD_M):}
#'       \deqn{ AD_M = n\sum_{x = 0}^{r} {\frac{{{{\left[ {{F_n}\left( x \right) - {F_{{{\hat \rho }_{X}(x)}}}\left( {x;\mu ,\phi } \right)} \right]}^2}{{\hat \rho }_{X}(x)}}}{{{F_{{{\hat \rho }_{X}(x)}}}\left( {x;\mu ,\phi } \right)\left[ {1 - {F_{{{\hat \rho }_{X}(x)}}}\left( {x;\mu ,\phi } \right)} \right]}}} }
#'     \item \strong{Probability Distance (PD):} Measures the sum of squared differences between
#'       empirical probabilities and the Stein-based estimated probabilities.
#'       \deqn{ PD = \sum_{j = 0}^{r} {\left( {{p_n}(j) - {{\hat \rho }_X}\left( j \right)} \right)} ^2}
#'   }
#'   Where \eqn{n} is the sample size, \eqn{r} is the maximum observed value, \eqn{p_n(x)} is the empirical probability,
#'   \eqn{F_n(x)} is the empirical distribution function, \eqn{p(x;\mu,\phi)} and \eqn{F(x;\mu,\phi)} are the
#'   mCMP PMF and CDF respectively, and \eqn{\hat{\rho}_X(x)} and \eqn{F_{\hat{\rho}_X}(x)} are the Stein-based
#'   estimated PMF and its corresponding CDF.
#'
#' @param sample Numeric sample data.
#' @param mu Numeric location parameter \eqn{\mu} (typically an estimated value).
#' @param phi Numeric dispersion parameter \eqn{\phi} (typically an estimated value).
#' @return A named numeric vector of test statistics (CV, AD, chi_sq, CV_M, AD_M, PD).
#' @seealso \code{\link{dcmp}} for the mCMP PMF, \code{\link{estimate_rho_mCMP}} for the Stein-based estimator.
#' @examples
#' # Simulate some mCMP data
#' set.seed(456)
#' sample_data <- simulate_mcmp_sample(n = 100, mu = 6, phi = 0.2)
#'
#' # Estimate parameters for the simulated data (as done in goodness_of_fit_mCMP)
#' estimated_params <- estimate_mu_phi(sample_data)
#' mu_est <- estimated_params[1]
#' phi_est <- estimated_params[2]
#'
#' # Compute all goodness-of-fit test statistics
#' stats <- compute_mod_statistics(sample = sample_data, mu = mu_est, phi = phi_est)
#' print(stats)
compute_mod_statistics <- function(sample, mu, phi) {
  n <- length(sample)
  r <- max(sample)
  freq_table <- table(factor(sample, levels = 0:r))
  Y_x <- as.numeric(freq_table)
  p_n <- Y_x / n
  F_n <- cumsum(p_n)

  hat_rho <- compute_hat_rho_vec(sample, mu, phi, r)
  F_rho <- cumsum(hat_rho)


  pmf <- dcmp(0:r, mu, phi)
  cdf <- cumsum(pmf)
  CV <- n * sum((F_n - cdf)^2 * pmf, na.rm = TRUE)

  CV_M <- n * sum((F_n - F_rho)^2 * hat_rho, na.rm = TRUE)

  denominator <- cdf * (1 - cdf)
  valid_idx <- which(denominator > 0)
  AD <- n * sum(((F_n[valid_idx] - cdf[valid_idx])^2 * pmf[valid_idx]) / denominator[valid_idx], na.rm = TRUE)



  denominator <- F_rho * (1 - F_rho)
  valid_idx <- which(denominator > 0)
  AD_M <- n * sum(((F_n[valid_idx] - F_rho[valid_idx])^2 * hat_rho[valid_idx]) / denominator[valid_idx], na.rm = TRUE)

  PD <- sum((p_n - hat_rho)^2)

  expected <- n * pmf
  chi_sq <- sum((Y_x - expected)^2 / expected, na.rm = TRUE)

  return(c(CV=CV, AD=AD, chi_sq = chi_sq, CV_M = CV_M, AD_M = AD_M, PD = PD))
}

#' @title Objective Function for Estimation (Internal Helper)
#'
#' @description This internal helper function defines the objective function \eqn{S(\mu,\phi)}
#'   that is minimized to estimate the mCMP parameters \eqn{\mu} and \eqn{\phi}.
#'   The estimation method used in the package is based on minimizing the probability distance (PD)
#'   between the empirical probabilities and the probabilities derived from Stein's characterization.
#'   The objective function is given by:
#'   \deqn{ S(\mu,\phi) = \sum_{j = 0}^{r} {\left( {p_n(j) - \hat{\rho}_X(j)} \right)} ^2 }
#'   where \eqn{p_n(j)} is the empirical probability of observing value \eqn{j} in the sample,
#'   and \eqn{\hat{\rho}_X(j)} is the estimated Stein-based probability for \eqn{j}
#'   (see \code{\link{estimate_rho_mCMP}}). This approach is computationally efficient as
#'   \eqn{\hat{\rho}_X(j)} does not involve the complex normalizing constant.
#'
#' @param par A numeric vector of two elements: \code{par[1]} is \eqn{\mu} and \code{par[2]} is \eqn{\phi}.
#' @param x A numeric vector of observed sample data.
#'
#' @return A numeric value representing the objective function value. Returns 0 if parameters lead to invalid calculations
#'   (e.g., \eqn{\mu \le 0} or issues with the adjustment term, as per the original code's logic).
#' @examples
#' # Simulate some mCMP data
#' set.seed(789)
#' sample_data <- simulate_mcmp_sample(n = 50, mu = 4, phi = -0.1)
#'
#' # Evaluate the objective function for some parameter guesses
#' # A guess close to the true parameters should yield a smaller value
#' S_val1 <- S_objective(par = c(mu = 4, phi = -0.1), x = sample_data)
#' cat("S_objective for (mu=4, phi=-0.1):", S_val1, "\n")
#'
#' # A guess further from the true parameters should yield a larger value
#' S_val2 <- S_objective(par = c(mu = 4.5, phi = 0), x = sample_data)
#' cat("S_objective for (mu=4.5, phi=0):", S_val2, "\n")
S_objective <- function(par, x) {
  mu <- par[1]
  phi <- par[2]
  check <- mu+((exp(phi) - 1) / (2 * exp(phi)))
  if (mu<=0||check <= 0|| is.na(check)) return(0)
  n <- length(x)
  r <- max(x)
  freq_table <- table(factor(x, levels = 0:r))
  p_n <- as.numeric(freq_table) / n
  hat_rho <- compute_hat_rho_vec(x, mu, phi, r)
  pn <- sum((p_n - hat_rho)^2)
  if (!is.finite(pn)){
    pn <- 0
  }
  return(pn)
}

#' @title Estimate mCMP Parameters by Minimizing S(mu, phi) (Internal Helper)
#'
#' @description This internal helper function estimates the parameters \eqn{\mu} and \eqn{\phi} of the
#'   mCMP distribution from observed data. It achieves this by numerically minimizing the
#'   \code{\link{S_objective}} function, which is based on the probability distance criterion.
#'   The minimization is performed using the "L-BFGS-B" method from the \code{\link[stats]{optim}} function.
#'   Initial guesses for \eqn{\mu} and \eqn{\phi} are derived from the sample mean and variance,
#'   respectively, to provide reasonable starting points for the optimization algorithm.
#'
#' @param x A numeric vector of observed sample data.
#' @return A numeric vector containing the estimated \eqn{\mu} and \eqn{\phi} values.
#' @examples
#' # Simulate some mCMP data
#' set.seed(101)
#' sample_data <- simulate_mcmp_sample(n = 70, mu = 7, phi = 0.5)
#'
#' # Estimate the parameters
#' estimated_params <- estimate_mu_phi(sample_data)
#' cat("Estimated mu:", estimated_params[1], "\n")
#' cat("Estimated phi:", estimated_params[2], "\n")
estimate_mu_phi <- function(x) {
  mu_start <- mean(x) # Initial guess for mu
  phi_start <- log(mean(x)/var(x)) # Initial guess for phi based on variance
  opt <- optim(
    par = c(mu_start, phi_start),
    fn = S_objective,
    x = x,
    method = "L-BFGS-B",
    lower = c(0, -Inf),
    upper = c(Inf, Inf),
    control = list(fnscale = 1)
  )
  return(opt$par)
}
