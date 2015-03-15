### Estimate the price of different European options on multiple stocks ########


### 0) Setup ###################################################################

## load packages
require(qrng)

## set seed for reproducibility
set.seed(271)


### 1) Functions ###############################################################

##' @title Payoff Function for Basket Options
##' @param K strike of the option
##' @param N notional of the option
##' @param S0 a d-vector with stocks' current prices
##' @param S a (n, d)-matrix with stocks' final prices
##' @param optiontype string "Call" or "Put"
##' @return n-vector values of the basket option's payoff function
##' @author Mathieu Cambou, Marius Hofert, Christiane Lemieux
##' @note This is the undiscounted value of the payoff function
##'       evaluated at the stock prices at maturity.
payoff <- function(K, N, S0, S, type = c("call", "put"),
                   method = c("basket", "best", "worst"))
{
  stopifnot(K >= 0, N >= 0, S0 >= 0, S >= 0,
            length(S0) == ncol(S))
  type <- match.arg(type)
  method <- match.arg(method)
  perf <- switch(method,
                 "basket" = {
                     rowMeans(t(t(S)/S0))
                 },
                 "best" = {
                     apply(t(t(S)/S0), 1, max)
                 },
                 "worst" = {
                     apply(t(t(S)/S0), 1, min)
                 },
                 stop("Wrong 'method'"))
  N * pmax(0, if(type=="call") perf - K else K - perf)
}

##TODO:
##' @title Payoff Function: Worst-Of Option
##' @return n-vector values of the worst-of option's payoff function

##' @title Payoff Function: Best-Of Option
##' @return n-vector values of the best-of option's payoff function


#TODO: the below function can be generalized to:
# - allow for correlations between the GBMs
# - allow for sampling a full path to use it for path-dependent options
##' @title Produces samples of a d dimensional Geometric Brownian motion (GBM)
##'        from samples of uniforms, given daily parameters and a time horizon
##'        (in number of days)
##' @param u (nxd)-matrix sample of uniforms
##' @param S0 d-vector with initial GBMs levels
##' @param sigma d-vector containing the volatilities of the GBMs
##' @param mu d-vector containing the drifts of the GBMs
##' @param T time horizon of the samples (in number of days)
##' @return (nxd)-matrix of GBM samples
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
generate_GBM <- function(u, S0, sigma, mu, T)
{
  stopifnot(length(sigma) == length(r),
            length(sigma) == length(S0),
            T >= 0, all(sigma >= 0), all(mu >= 0))

  log.drift <- (mu - 0.5 * sigma * sigma)*T # (1xd)
  # log.diffusion <- sweep(qnorm(u), 2, sigma, '*') * sqrt(T) # (nxd)
  log.diffusion <- t(t(qnorm(u))*sigma)  # (nxd)
  log. <- log.drift + log.diffusion  # (nxd)
  t(t(exp(log.))*S0) # (nxd)
}


### 2) Case Study ##############################################################

## Simulation parameters
n <- 10^6 # Monte Carlo sample size
d <- 4 # dimension

## Stochastic process parameters
sigma <- rep(0.2, d) # volatilities
r. <- 0.0001 # continuously compounded short rate
r <- rep(r., d)
S0 <- rep(100, d) # initial stocks' levels
K <- 1.1 # option strike
N <- 1000 # option notional

## Generates Uniforms using Quasi
generator <- 17
u.quasi <- korobov(n, d = d, generator = generator) # TODO randomize
## => possibly copula-transform

## Generates GBM from Quasi Samples
S.quasi <- generate_GBM(u.quasi, S0, sigma, r, T)

## Call options using the Quasi samples
price.basketCall  <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Call"))
price.worstOfCall <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Call"), method="worst")
price.bestOfCall  <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Call"), method="best")

## Put options using the Quasi samples
price.basketPut  <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Put"))
price.worstOfPut <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Put"), method="worst")
price.bestOfPut  <- exp(-r.*T) * mean(payoff(K, N, S0, S.quasi, type = "Put"), method="best")

## Clayton (CDM + MO) / Gumbel (MO) / FGM (both)



## TODO: mathieu: insurance example