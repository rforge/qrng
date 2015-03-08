### Interfaces to quasi-random sequences #######################################


##' @title Korobov's quasi-random sequence
##' @param n number of points (>= 2 as generator has to be in {1,..,n-1}
##' @param d dimension
##' @param generator generator in {1,..,n-1}; either a vector of length d
##'        or a single number (which is appropriately extended)
##' @return an (n, d)-matrix containing the quasi-random sequence
##' @author Marius Hofert
korobov <- function(n, d, generator) {
    stopifnot(n >= 2, d >= 1, (l <- length(generator)) == 1 || l == d,
              1 <= generator, generator <= n-1, generator %% 1 == 0)
    if(l == 1) generator <- generator^(0:(d-1)) %% n # vectorize
    u <- .Call(korobov_, n, d, generator)
    if(d == 1) as.vector(u) else u
}

