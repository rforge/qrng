### Auxiliaries ################################################################

##' @title Convert Matrices to Arrays
##' @param x (n, f * d)-matrix, e.g. sobol(n, d = f * d, randomize = "digital.shift")
##' @param f integer f >= 1 dividing ncol(x)
##' @param method method how output is formatted
##' @return an array of the form as given by 'method'
##' @author Marius Hofert
##' @note This is helpful for time series applications where f = number of time
##'       steps to simulate
to_array <- function(x, f, method = c("(n*f,d)", "(n,f,d)"))
{
    stopifnot(is.matrix(x), f >= 1)
    if(f == 1) return(x)
    dm <- dim(x)
    n  <- dm[1]
    d. <- dm[2] # = f * d
    if(d. %% f != 0)
        stop("'f' must divide d = ncol(x)")
    d <- d./f
    method <- match.arg(method)
    switch(method,
           "(n*f,d)" = { # convert (n, f * d)-matrix into (n * f, d)-matrix
               matrix(t(x), ncol = d, byrow = TRUE)
           },
           "(n,f,d)" = { # convert (n, f * d)-matrix into (n, f, d)-array
               array(x, dim = c(n, f, d))
           },
           stop("Wrong 'method'"))
}
