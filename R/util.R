# Random variate generation from Laplace distribution
rlapl <- function(n, mu, lambda) {
  u <- runif(n)
  cond <- (u >= 0.5)
  rand <- mu + (-1)^cond * lambda * log(2 * (cond + (-1)^cond * u))
  return(rand)
}

# Probability density function (pdf) of Laplace distribution
dlapl <- function(x, mu = 0, lambda = 1) {
  pdf <- 1 / (2 * lambda) * exp(-abs(x - mu) / lambda)
  return(pdf)
}

# Cumulative distribution function (cdf) of Laplace distribution
plapl = function(x, mu = 0, lambda = 1) {
  cond <- (x >= mu)
  cdf <- cond + (-1)^cond*0.5 * exp(-abs(x - mu) / lambda)
  return(cdf)
}

# Log-likelihood (normal error)
#' @importFrom stats pnorm
llik.n <- function(a, s, x) {
  n <- length(x)
  return(-n * log(2 * a) + sum(log(pnorm((x + a) / s) - pnorm((x - a) / s))))
}

#' @importFrom stats pnorm
llik.n2 = function(t, x) {
  n <- length(x)
  a <- t[1]
  s <- t[2]
  return (-n * log(2 * a) + sum(log(pnorm((x + a) / s)-pnorm((x - a) / s))))
}

# The function whose zeros are endpoints of LR confidence interval (normal
# error)
#' @importFrom stats qchisq
fun.cin <- function(a, ml, x, s, cl) {
  -llik.n(a, s, x) + llik.n(ml, s, x) - 0.5 * qchisq(cl, 1)
}

# The integrand appearing in the Fisher information (normal error)
#' @importFrom stats dnorm
#' @importFrom stats pnorm
int.n <- function(x, alfa) {
  (dnorm(x + alfa) + dnorm(x - alfa))^2 / (pnorm(x - alfa, lower.tail = F) -
                                           pnorm(x + alfa, lower.tail = F))
}

# Log-likelihood (Laplace error)
llik.l = function(a, l, x) {
  n <- length(x)
  return(-n * log(2 * a) + sum(log(plapl((x + a) / l) - plapl((x - a) / l))))
}

llik.l2 = function(t, x) {
  n <- length(x)
  a <- t[1]
  l <- t[2]
  return(-n * log(2 * a) + sum(log(plapl((x + a) / l) - plapl((x - a) / l))))
}

# The function whose zeros are endpoints of LR confidence interval (Laplace
# error)
#' @importFrom stats qchisq
fun.cil <- function(a, ml, x, l, cl) {
  -llik.l(a, l, x) + llik.l(ml, l, x) - 0.5 * qchisq(cl, 1)
}

# The integrand appearing in the Fisher information (Laplace error)
int.l <- function(x, alfa) {
  (dlapl(x + alfa) + dlapl(x - alfa))^2/(plapl(x + alfa) - plapl(x - alfa))
}

# Log-likelihood (scaled Student error with df = 5)
#' @importFrom stats pt
llik.s = function(a, s, x) {
  n <- length(x)
  return(-n * log(2 * a) + sum(log(pt((x + a) / s, 5) - pt((x - a) / s, 5))))
}

#' @importFrom stats pt
llik.s2 = function(t, x) {
  n <- length(x)
  a <- t[1]
  s <- t[2]
  return(-n * log(2 * a) + sum(log(pt((x + a) / s, 5) - pt((x - a) / s, 5))))
}

# The function whose zeros are endpoints of LR confidence interval (Student
# error)
#' @importFrom stats qchisq
fun.cis <- function(a, ml, x, s, cl) {
  -llik.s(a, s, x) + llik.s(ml, s, x) - 0.5 * qchisq(cl, 1)
}

# The integrand appearing in the Fisher information (Student error)
#' @importFrom stats dt
#' @importFrom stats pt
int.s <- function(x, alfa) {
  (dt(x + alfa, 5) + dt(x - alfa, 5))^2/(pt(x + alfa, 5) - pt(x - alfa, 5))
}
