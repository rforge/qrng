### Estimate the price of different European options on multiple stocks ########


### 0) Setup ###################################################################

## load packages
require(qrng)


### 1) Functions ###############################################################

##' @title Payoff Function for Basket Options
##' @param K strike of the option
##' @param N notional of the option
##' @param S0 a d-vector with stocks' current prices
##' @param S a (n, d)-matrix with stocks' final prices
##' @param type string "call" or "put"
##' @param method option type; available are "basket", "worst-of" and "best-of"
##'        option
##' @return n-vector values of the basket option's payoff function
##' @author Mathieu Cambou, Marius Hofert, Christiane Lemieux
##' @note This is the undiscounted value of the payoff function
##'       evaluated at the stock prices at maturity.
payoff <- function(K, N, S0, S, type = c("call", "put"),
                   method = c("basket", "worst.of", "best.of"))
{
  stopifnot(K >= 0, N >= 0, S0 >= 0, S >= 0,
            length(S0) == ncol(S))
  type <- match.arg(type)
  method <- match.arg(method)
  perf <- switch(method,
                 "basket" = {
                     rowMeans(t(t(S)/S0))
                 },
                 "worst.of" = {
                     apply(t(t(S)/S0), 1, min)
                 },
                 "best.of" = {
                     apply(t(t(S)/S0), 1, max)
                 },
                 stop("Wrong 'method'"))
  N * pmax(0, if(type=="call") perf - K else K - perf)
}

##' @title Produces Samples of a d-dimensional Geometric Brownian Motion (GBM)
##' @param u (n,d)-matrix sample of uniforms
##' @param S0 d-vector with initial GBM levels
##' @param sigma d-vector containing the (daily) volatilities of the GBMs
##' @param mu d-vector containing the drifts of the GBMs
##' @param T time horizon of the samples (in number of days)
##' @return (n,d)-matrix of GBM samples (d independent BMs)
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
##' @note Further generalizations are possible:
##'       - allow for correlations between the GBMs
##'       - sample a full path for pricing path-dependent options
rGeoBM <- function(u, S0, sigma, mu, T)
{
  stopifnot(0 < u, u < 1, (d <- length(S0)) == length(sigma), sigma >= 0,
            length(mu) == d, mu >= 0, r >= 0, T >= 0)
  log.diffusion <- qnorm(u) * matrix(rep(sigma, n), ncol=d, byrow=TRUE) # (n,d)-matrix; or t(t(qnorm(u))*sigma)
  log.drift <- (mu - sigma^2 / 2) * T # d-vector
  log. <- matrix(rep(log.drift, n), ncol=d, byrow=TRUE) + log.diffusion # (n,d)-matrix
  matrix(rep(S0, n), ncol=d, byrow=TRUE) * exp(log.) # S_t, t in 1,..,T; (n,d)-matrix
}


### 2) Case Study ##############################################################

## Simulation parameters
n <- 1e5 # Monte Carlo sample size
d <- 4 # dimension

## Stochastic process parameters
sigma <- rep(0.2, d) # volatilities
r. <- 0.0001 # continuously compounded short rate
r <- rep(r., d)
S0 <- rep(100, d) # initial stocks' levels
K <- 1.1 # option strike
N <- 1000 # option notional

## Generates generalized Halton sequence
set.seed(271)
u.quasi <- ghalton(n, d = d)
## => possibly copula-transform

## Generates Geometric Brownian Motion from quasi-random sequence
S.quasi <- rGeoBM(u.quasi, S0, sigma, r, T)

## Call options using the quasi-random sequence
basket.call  <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "call"))
worst.of.call <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "call", method="worst.of"))
## TODO: the same?

## TODO: Mathieu, did you have an insurance example in mind? maybe a simple one which runs
##       fast enough to be part the demo (like around max. 30s)

## TODO: need to write a script (not this demo) for the following simulation study:
## - QRNG: ghalton, sobol, none
## - dimensions d: 5, 20
## - copulas: Clayton, t
## - construction methods: CDM, MO
## - Kendall's tau: 0.2, 0.5
## - bootstrap replications: 25 (different randomizations => estimate variance)
## - options: 'basket' and 'worst-of' (or 'best-of'?)
## - a bunch of different 'n': 64 different n (for creating plots in 'n')
##   => generate for the largest sample size and then take sub-parts
