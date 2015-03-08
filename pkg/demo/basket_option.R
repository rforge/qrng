### Estimate the price of a basket put option ##################################

require(qrng)

n <- 2
d <- 1
generator <- 1
(u <- korobov(n, d=d, generator=generator))


pairs(u, gap=0, pch=".", labels=as.expression(
      sapply(1:d, function(j) bquote(italic(u[.(j)])))))

