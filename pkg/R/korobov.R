### Main interface to (randomized) quasi-random number generators ##############

# TODO write man

##' @title Korobov's quasi-random sequence
##' @param n number of points
##' @param d dimension
##' @param generator generator; either a vector of length d or a single number
##' @return an (n, d)-matrix containing the quasi-random sequence
##' @author Marius Hofert
##' @note - n should be a prime or a power of two
##'       - if n is a power of two, the generator has to be odd
##'       - the generator should be a primitive element mod n ()
korobov <- function(n, d, generator) {
    stopifnot(n >= 1, d >= 1, 1 <= generator, generator <= n-1)
    l <- length(generator)
    if(l != 1 && l != d)
        stop("'generator' must be of length 1 or ", d)
    if(l == 1) {
        stopifnot(generator %% 1 == 0) # check if integer
        generator <- generator^(0:(d-1)) %% n
    } else stopifnot(generator %% 1 == 0)
    .Call(korobov_, n, d, generator)
}
