### Interfaces to quasi-random sequences #######################################


##' @title Korobov's sequence
##' @param n number of points (>= 2 as generator has to be in {1,..,n-1}
##' @param d dimension
##' @param generator generator in {1,..,n-1}; either a vector of length d
##'        or a single number (which is appropriately extended)
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
korobov <- function(n, d, generator)
{
    stopifnot(n >= 2, d >= 1, (l <- length(generator)) == 1 || l == d,
              1 <= generator, generator <= n-1, generator %% 1 == 0)
    lim <- 2^31-1
    if(n > lim)
        stop("'n' has to be <= 2^31-1")
    if(d > lim)
        stop("'d' has to be <= 2^31-1")
    if(l == 1) generator <- generator^(0:(d-1)) %% n # vectorize
    u <- .Call(korobov_, n, d, generator)
    if(d == 1) as.vector(u) else u
}

##' @title Generalized Halton sequence
##' @param n number of points
##' @param d dimension
##' @param n.random
##' @return an (n, d)-matrix (an n-vector if d=1) containing the
##'         quasi-random sequence
##' @author Marius Hofert
ghalton <- function(n, d, n.random)
{
    stopifnot(n >= 1, d >= 1, n.random >= 0)
    if(n > 2^32-1)
        stop("'n' has to be <= 2^32-1")
    if(d > 360)
        stop("'d' has to be <= 360")
    if(n.random > 2^31-1)
        stop("'n.random' has to be <= 2^31-1")
    ghalton_ <- NULL # to make CRAN check happy (for some reason not required for korobov())
    u <- .Call(ghalton_, n, d, n.random)
    if(d == 1) as.vector(u) else u
}