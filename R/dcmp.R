#' @title The Conway-Maxwell Poisson Probability Mass Function
#'
#' @description Computes the probability mass function (PMF) of the Conway-Maxwell Poisson (mCMP) distribution.
#'   The PMF is given by:
#'   \deqn{ p(X=x;\mu,\phi) = \frac{\left( \mu + \frac{e^\phi - 1}{2e^\phi} \right)^{x e^\phi}}{(x!)^{e^\phi} Z(\mu,\phi)} }
#'   where \eqn{x = 0, 1, 2, \ldots}, \eqn{\mu > 0}, \eqn{-\infty < \phi < \infty}, and \eqn{Z(\mu,\phi)} is the normalizing constant.
#'   The normalizing constant is defined as:
#'   \deqn{ Z(\mu,\phi) = \sum_{j=0}^{\infty} \frac{\left( \mu + \frac{e^\phi - 1}{2e^\phi} \right)^{j e^\phi}}{(j!)^{e^\phi}} }
#'   This function uses an iterative summation for \eqn{Z(\mu,\phi)} up to a specified number of \code{max_terms}
#'   or until the added terms fall below a given \code{tol}.
#'
#' @param x A vector of non-negative integers.
#' @param mu A numeric value representing the location parameter \eqn{\mu}. Must be positive.
#' @param phi A numeric value representing the dispersion parameter \eqn{\phi}.
#' @param max_terms An integer for the maximum number of terms to sum for the normalizing constant.
#' @param tol A numeric tolerance for the normalizing constant sum. Summation stops when
#'   the ratio of the current term to the total sum falls below this tolerance.
#'
#' @return A numeric vector of probabilities.
#' @examples
#' dcmp(0:5, mu = 2, phi = 0) # Equi-dispersed (Poisson)
#' dcmp(0:5, mu = 2, phi = 1) # Under-dispersed
#' dcmp(0:5, mu = 2, phi = -1) # Over-dispersed
dcmp <- function(x, mu, phi, max_terms = 1000, tol = 1e-12) {
  # Input validation
  if (any(x < 0) || any(x != floor(x))) stop("x must be non-negative integers")
  if (mu < 0) stop("mu must be positive")

  # Normalizing constant
  Z_val <- compute_Z(mu, phi, max_terms, tol)
  if (!is.finite(Z_val) || Z_val <= 0) return(rep(NA_real_, length(x)))

  # Adjustment term
  adjustment <- (exp(phi) - 1) / (2 * exp(phi))
  base_term <- mu + adjustment
  if (!is.finite(base_term) || base_term <= 0) return(rep(0, length(x)))

  result <- numeric(length(x))

  for (i in seq_along(x)) {
    if (x[i] == 0) {
      log_numerator <- 0
      log_denominator <- 0
    } else {
      log_numerator <- x[i] * exp(phi) * log(base_term)
      log_denominator <- exp(phi) * lgamma(x[i] + 1) # <- stable factorial
    }

    # log(PMF) = log_num - log_denom - log(Z)
    log_p <- log_numerator - log_denominator - log(Z_val)
    result[i] <- exp(log_p)
  }

  return(result)
}

#' @title Compute Normalizing Constant (Internal)
#' @description Internal helper function to compute the normalizing constant Z for the mCMP distribution.
#' @param mu Numeric location parameter.
#' @param phi Numeric dispersion parameter.
#' @param max_terms Maximum terms to sum.
#' @param tol Tolerance for series convergence.
#' @return Numeric value of Z.
compute_Z <- function(mu, phi, max_terms = 1000, tol = 1e-12) {
  adjustment <- (exp(phi) - 1) / (2 * exp(phi))
  base_term <- mu + adjustment

  if (!is.finite(base_term) || base_term <= 0) return(NaN)

  Z_sum <- 0
  for (j in 0:max_terms) {
    if (j == 0) {
      term <- 1
    } else {
      log_numerator <- j * exp(phi) * log(base_term)
      log_denominator <- exp(phi) * lgamma(j + 1) # <- stable factorial

      if (log_numerator - log_denominator < log(tol)) break

      term <- exp(log_numerator - log_denominator)
    }

    if (!is.finite(term)) break
    Z_sum <- Z_sum + term

    if (j > 0 && is.finite(term / Z_sum) && term / Z_sum < tol) break
  }

  return(Z_sum)
}
