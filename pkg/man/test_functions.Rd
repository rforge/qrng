\name{test_functions}
\alias{sum_of_squares}
\alias{sobol_g}
\alias{exceedance}
\title{Test Functions}
\description{
  Functions for testing low-discrepancy sequences.
}
\usage{
sum_of_squares(u)
sobol_g(u, copula = indepCopula(dim = ncol(u)), alpha = 1:ncol(u), ...)
exceedance(x, p = 0.99, method = c("indicator", "individual.given.sum.exceeds",
                                   "sum.given.sum.exceeds"))
}
\arguments{
  \item{u}{\eqn{(n, d)}-matrix containing \eqn{n}
    \eqn{d}-dimensional realizations (of a potential quasi-random number
    generator). For \code{sum_of_squares()} these need to be marginally
    standard uniform and for \code{sobol_g()} they need to follow the
    copula specified by \code{copula}.}
  \item{copula}{\code{\link{Copula}} object for which the inverse
    Rosenblatt transformation exists.}
  \item{alpha}{vector of parameters of Sobol's g test function.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{cCopula}()}.}
  \item{x}{\eqn{(n, d)}-matrix containing \eqn{n} \eqn{d}-dimensional
    realizations.}
  \item{p}{
    \describe{
      \item{"indicator"}{\eqn{d}-vector containing the probabilities
	determining componentwise thresholds via empirical quantiles.}
      \item{"individual.given.sum.exceeds", "sum.given.sum.exceeds"}{
	probability determining the threshold for the sum (row sums of
	\code{x}) via the corresponding empirical quantile.}
    }
  }
  \item{method}{\code{\link{character}} string indicating the type of
    exceedance computed (see Section Value below).}
}
\value{
  \code{sum_of_squares()} returns an \eqn{n}-vector
  (\code{\link{numeric}(n)}) with the rowwise computed scaled sum
  of squares (theoretically integrating to 1).

  \code{sobol_g()} returns an \eqn{n}-vector (\code{\link{numeric}(n)})
  with the rowwise computed Sobol' g functions.

  \code{exceedance()}'s return value depends on \code{method}:
  \describe{
    \item{"indicator"}{returns indicators whether,
      componentwise, \code{x} exceeds the threshold determined as the
      \code{p}-quantile.}
    \item{"individual.given.sum.exceeds"}{returns all rows of \code{x}
      whose sum exceeds the threshold determined as the
      \code{p}-quantile of the row sums of \code{x}.}
    \item{"sum.given.sum.exceeds"}{returns the row sums of those
      rows of \code{x} whose sum exceeds the threshold determined as the
      \code{p}-quantile of the row sums of \code{x}.}
  }
}
\details{
  See \code{VaR_np(<matrix>)} and \code{ES_np(<matrix>)} from \pkg{qrmtools}
  for two more test functions.
}
\author{Marius Hofert and Christiane Lemieux}
\references{
  Radovic, I., Sobol, I. M., and Tichy, R. F. (1996). Quasi-Monte Carlo
  methods for numerical integration: Comparison of different low
  discrepancy sequences. \emph{Monte Carlo Methods and Applications}
  \bold{2}(1), 1--14.

  Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
  A Comparative Study. \emph{ACM-TOMACS} \bold{19}(4), Article 15.
}
\examples{
n <- 1021
d <- 4
set.seed(271)
u <- sobol(n, d, randomize = "Owen") # Sobol' randomized with Owen scrambling
mean(sum_of_squares(u)) # estimate of E(3(sum_{j=1}^d U_j^2)/d)
mean(sobol_g(u)) # estimate of E(<Sobol's g function>)
mean(exceedance(qnorm(u))) # estimate of E(<exceedance indicator over q>) = P(X > q)
}
\keyword{utilities}