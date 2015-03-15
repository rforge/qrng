### Experimenting with test functions for MC ###################################

## NOTE: this is an experimental playground...

## see Sobol, Asotsky (2003, "One more experiment on estimating high-dimensional
## integrals by quasi-Monte Carlo methods"), Mathematics and Computers in
## Simulation, 62 (3--6), 255--263

require(qrng)
require(copula)
require(randtoolbox)

## Test function
sobolAsot <- function(x, d, a) {
    apply(sapply(1:d, function(j) (abs(4*x[,j]-2) + a[j])/(1+a[j])), 1, prod)
}

## Korobov's sequence
n <- 1021 # number of points
d <- 5 # dimension
generator <- 76 # generator
U1 <- korobov(n, d=d, generator=generator) # Korobov's sequence

## Halton's sequence
U2 <- halton(n, dim=d)

## Copula sample
family <- "Gumbel" # copula family
tau <- 0.5 # Kendall's tau
th <- iTau(getAcop(family), tau) # corresponding copula parameter
cop <- onacopulaL(family, nacList=list(th, 1:d)) # define copula
set.seed(271) # set seed
U3 <- rCopula(1e5, cop) # sample from cop

## Evaluate test function
U <- U1
mean(sobolAsot(U, d=d, a=rep(0.01, d)))
mean(sobolAsot(U, d=d, a=rep(1, d)))
mean(sobolAsot(U, d=d, a=1:d))
mean(sobolAsot(U, d=d, a=(1:d)^2))
mean(sobolAsot(U, d=d, a=(d- 1:d + 1)^2))
