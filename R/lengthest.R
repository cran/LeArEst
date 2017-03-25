#' Performs width estimation for a numerical data set.
#'
#' Function \code{lengthest()} computes the length of an interval which is the
#' domain of uniform distribution from data contaminated with additive error.
#' The additive error can be chosen as Laplace, Gauss or scaled Student
#' distribution with 5 degrees of freedom.
#'
#' @param x Vector of input data.
#' @param error A character string specifying the error distribution. Must be
#'   one of "laplace", "gauss" or "student". Can be abbreviated.
#' @param var Explicit error variance. Needs to be given if \code{var.est} is
#'   not given.
#' @param var.est A character string specifying the method of error variance
#'   estimation. Must be given if \code{var} is not given. Can be "MM" (Method
#'   of Moments) or "ML" (Maximum Likelihood).
#' @param conf.level Confidence level of the confidence interval.
#'
#' @return List containing:
#'   \itemize{
#'     \item error.type: A character string describing the type of the
#'       error distribution.
#'     \item radius: Estimated half-width of uniform distribution.
#'     \item var.error: Error variance, estimated or given.
#'     \item conf.level: Confidence level of the confidence interval.
#'     \item method: A character string indicating what method for
#'       computing a confidence interval was used (asymptotic distribution of
#'       ML or likelihood ratio statistic).
#'     \item conf.int: The confidence interval for half-width.
#'   }
#'
#' @examples
#' # generate uniform data with additive error and run a length estimation on it
#' sample_1 <- runif(1000, -1, 1)
#' sample_2 <- rnorm(1000, 0, 0.1)
#' sample <- sample_1 + sample_2
#' out <- lengthest(x = sample, error = "gauss", var.est = "MM", conf.level = 0.90)
#'
#' @importFrom stats sd
#' @importFrom stats optim
#' @importFrom stats integrate
#' @importFrom stats qnorm
#' @importFrom stats uniroot
#'
#' @export
lengthest <- function(x,
                      error = c("laplace", "gauss", "student"),
                      var = NULL,
                      var.est = c("MM", "ML"),  #-- Safet origigi
                      conf.level = 0.95) {
  cl <- conf.level
  c1 <- list(fnscale = -1)
  m2 <- mean(x^2)
  m4 <- mean(x^4)
  error <- match.arg(error)

  n <- length(x)

  lower <- 0
  start <- max(abs(x))
  upper <- 2*start

  options(warn = -1)

  if (error == "laplace") {
    if (is.null(var) & var.est == "ML") {
      start = c(max(abs(x)),sd(x))
      as.ml = optim(start, llik.l2, x = x, control = c1)[[1]]
      a.ml <- as.ml[1]
      l <- as.ml[2]
    } else {
      if (is.null(var) & var.est == "MM") {
        l <- 1/sqrt(6)*sqrt(-2*m2 + sqrt(5*(m4 - m2^2)))
        # l<-sqrt(-1)
        if (is.na(l)) {
          stop("MM estimate of error variance doesn't exist.
  Try again using ML estimate or specific value for error variance.")
          }
        } else {
          l <- sqrt(var/2)
        }
        a.ml <- optim(start, llik.l, x = x, l = l, method = "Brent",
          lower = lower, upper = upper, control = c1)[[1]]
    }
    alfa <- a.ml/l
    if (fun.cil(1e-06,a.ml,x,l,cl)*
        fun.cil(a.ml,a.ml,x,l,cl)>0 ||
        fun.cil(a.ml,a.ml,x,l,cl)*
        fun.cil(1e+06,a.ml,x,l,cl)>0) {
      method <- "ADMLE"
      integral <- integrate(int.l, alfa = alfa, 0, alfa)[[1]] +
                  (cosh(alfa))^2/sinh(alfa)*exp(-alfa)
      avar.ml <- 1/(n*(-1/a.ml^2+1/(a.ml*l)*integral))
      r <- qnorm(0.5*(1 + cl))*sqrt(avar.ml)
      ci <- c(a.ml - r, a.ml + r)
    } else {
      method <- "ADLR"
      l.int <- uniroot(fun.cil, ml = a.ml, x = x, l = l, cl = cl,
                       lower = 1e-06, upper = a.ml)$root
      r.int <- uniroot(fun.cil, ml = a.ml, x = x, l = l, cl = cl,
                       lower = a.ml, upper = 1e+06)$root
      ci <- c(l.int,r.int)
    }
    v <- 2*l^2
    } else if (error == "gauss") {
      if (is.null(var) & var.est == "ML") {
        start = c(max(abs(x)),sd(x))
        as.ml = optim(start, llik.n2, x = x, control = c1)[[1]]
        a.ml <- as.ml[1]
        s <- as.ml[2]
      } else {
        if (is.null(var) & var.est == "MM") {
          s <- sqrt(mean(x^2)-sqrt(5/6*(3*(mean(x^2))^2-mean(x^4))))
          if (is.na(s)) {
            stop("MM estimate of error variance doesn't exist.
  Try again using ML estimate or specific value for error variance.")
          }
          } else {
            s <- sqrt(var)
          }
          a.ml <- optim(start, llik.n, x = x, s = s, method = "Brent",
                        lower = lower, upper = upper, control = c1)[[1]]
        }
      alfa <- a.ml/s
      if (fun.cin(1e-06,a.ml,x,s,cl)*
          fun.cin(a.ml,a.ml,x,s,cl)>0 ||
          fun.cin(a.ml,a.ml,x,s,cl)*
          fun.cin(1e+06,a.ml,x,s,cl)>0) {
        method <- "ADMLE"
        integral <- integrate(int.n, alfa = alfa , 0, alfa  + 37.5)[[1]]
        avar.ml <- 1/(n*(-1/a.ml^2 + 1/(a.ml*s)*integral))
        r <- qnorm(0.5*(1 + cl))*sqrt(avar.ml)
        ci <- c(a.ml - r, a.ml + r)
      } else {
        method <- "ADLR"
        l.int <- uniroot(fun.cin, ml = a.ml, x = x, s = s, cl = cl,
                         lower = 1e-06, upper = a.ml)$root
        r.int <- uniroot(fun.cin, ml = a.ml, x = x, s = s, cl = cl,
                         lower = a.ml, upper = 1e+06)$root
        ci <- c(l.int,r.int)
      }
      v <- s^2
      } else {
        if (is.null(var) & var.est == "ML") {
          start = c(max(abs(x)),sd(x))
          as.ml = optim(start, llik.s2, x = x, control = c1)[[1]]
          a.ml <- as.ml[1]
          s <- as.ml[2]
        } else {
          if (is.null(var) & var.est == "MM") {
            s <- 1/(2*sqrt(5))*sqrt(-3*m2 + sqrt(15*(2*m4 - 3*m2^2)))
            if (is.na(s)) {
              stop("MM estimate of error variance doesn't exist.
  Try again using ML estimate or specific value for error variance.")
              }
            } else {
              s <- sqrt(var*3/5)
            }
            a.ml <- optim(start, llik.s, x = x, s = s, method = "Brent",
                          lower = lower, upper = upper, control = c1)[[1]]
        }
        alfa <- a.ml/s
        if (fun.cis(1e-06,a.ml,x,s,cl)*
            fun.cis(a.ml,a.ml,x,s,cl)>0 ||
            fun.cis(a.ml,a.ml,x,s,cl)*
            fun.cis(1e+06,a.ml,x,s,cl)>0) {
          method <- "ADMLE"
          integral <- integrate(int.s, alfa = alfa, 0, alfa + 200)[[1]]
          avar.ml <- 1/(n*(-1/a.ml^2+1/(a.ml*s)*integral))
          r <- qnorm(0.5*(1 + cl))*sqrt(avar.ml)
          ci <- c(a.ml - r, a.ml + r)
        } else {
          method <- "ADLR"
          l.int <- uniroot(fun.cis, ml = a.ml, x = x, s = s, cl = cl,
                           lower = 1e-06, upper = a.ml)$root
          r.int <- uniroot(fun.cis, ml = a.ml, x = x, s = s, cl = cl,
                           lower = a.ml, upper = 1e+06)$root
          ci <- c(l.int,r.int)
        }
        v <- s^2*5/3
        }
  if (!is.null(var)) {
    names(v) <- "known error variance:"
  } else if (var.est == "ML") {
    names(v) <- "ML estimate for error variance:"
  } else {
    names(v) <- "MM estimate for error variance:"
  }
  if (method == "ADLR") {
    method <- "Asymptotic distribution of LR statistic"
  } else {
    method <- "Asymptotic distribution of MLE"
  }
  names(a.ml) <- "MLE for radius (a) of uniform distr.:"

  options(warn = 0)

  out <- list(error.type = error,
    radius = a.ml,
    var.error = v,
    conf.level = cl,
    method = method,
    conf.int = ci)

  return(out)
}
