### Estimate the price of a basket put option

require(qrng)

n <- 2^10
d <- 4
generator <- 17
u <- korobov(n, d=d, generator=generator)
pairs(u, gap=0)

## TODO:
## - n should be a prime or a power of two
## - if n is a power of two, the generator has to be odd
## - the generator should be a primitive element mod n ()
## - send her standalone file for korobov => need that for gHalton