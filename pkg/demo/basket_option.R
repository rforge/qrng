### Estimate the price of a basket put option ##################################

require(qrng)

n <- 2^10
d <- 4
generator <- 17
u <- korobov(n, d=d, generator=generator)
pairs(u, gap=0, pch=".", labels=as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

