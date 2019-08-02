## By Marius Hofert

## Examples of test function calls

## Generate some (here: copula, pseudo-random) data
library(copula)
set.seed(271)
cop <- claytonCopula(iTau(claytonCopula(), tau = 0.5)) # Clayton copula
U <- rCopula(1000, copula = cop)

## Compute sum of squares test function
mean(sum_of_squares(U)) # estimate of E(3(sum_{j=1}^d U_j^2)/d)

## Compute the Sobol' g test function
if(packageVersion("copula") >= "0.999-20")
    mean(sobol_g(U)) # estimate of E(<Sobol's g function>)

## Compute an exceedance probability
X <- qnorm(U)
mean(exceedance(X, q = qnorm(0.99))) # fixed threshold q
mean(exceedance(X, p = 0.99)) # empirically estimated marginal p-quantiles as thresholds

## Compute 99% expected shortfall for the sum
mean(exceedance(X, p = 0.99, method = "sum.given.sum.exceeds")) # version 1
library(qrmtools)
ES_np(X, level = 0.99) # version 2
