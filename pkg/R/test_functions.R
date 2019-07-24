### Test functions #############################################################

##' @title Sum of Squares Test Function
##' @param u (n, d)-matrix containing the n d-dimensional realizations marginally
##'        standard uniform
##' @return standardized sum of squares function (integrates to 1)
##' @author Marius Hofert
##' @note If you want an estimator of E(<test function>) call mean(sum_of_squares())
sum_of_squares <- function(u)
{
    if(!is.matrix(u)) u <- rbind(u)
    3 * rowSums(u^2) / ncol(u)
}

##' @title Sobol's g Test Function
##' @param u (n, d)-matrix containing the n d-dimensional realizations from
##'        copula 'copula'
##' @param copula a copula object
##' @param alpha parameters alpha of Sobol's g
##' @param ... additional arguments passed to the underlying cCopula()
##' @return Sobol's g
##' @author Marius Hofert
##' @note - See Radovic, Sobol, Tichy (1996, "Quasi-Monte Carlo Methods for
##'         Numerical Integration: Comparison of Different Low Discrepancy
##'         Sequences"
##'       - If you want an estimator of E(<test function>) call mean(sobol_g())
sobol_g <- function(u, copula = indepCopula(dim = ncol(u)), alpha = 1:ncol(u), ...)
{
    if(packageVersion("copula") < "0.999-20")
        stop('Your version of \'copula\' is not sufficient. Consider updating via install.packages("copula", repos = "http://R-Forge.R-project.org")')
    v <- cCopula(u, copula = copula, inverse = TRUE, ...)
    a <- rep(alpha, each = nrow(v))
    apply((abs(4 * v - 2) + a) / (1 + a), 1, prod)
}

##' @title Computing Indicators of Rows Exceeding a Threshold
##' @param x (n, d)-matrix containing the n d-dimensional realizations
##' @param p d-vector of probabilities determining the (joint) threshold
##'        through componentwise quantiles
##' @param method character string indicating the type of exceedance computed
##' @return Exceedance indicators whether X > x
##' @author Marius Hofert
##' @note - Run time can be improved for large d by checking the first dimension
##'         that falls below p.
##'       - If you want an estimator of P(X > x) = E(<test function>),
##'         call mean(exceedance())
exceedance <- function(x, p = 0.99, method = c("indicator",
                                               "individual.given.sum.exceeds",
                                               "sum.given.sum.exceeds"))
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x)
    method <- match.arg(method)
    switch(method,
           "indicator" = {
               q <- if(length(p) == 1) {
                        apply(x, 2, quantile, probs = p, names = FALSE, type = 1)
                    } else {
                        stopifnot(length(p) == d)
                        sapply(seq_len(d), function(j)
                            quantile(x[,j], probs = p[j], names = FALSE, type = 1))
                    }
               rowSums(x > rep(q, each = nrow(x))) == d
           },
           "individual.given.sum.exceeds" = {
               stopifnot(length(p) == 1)
               s <- rowSums(x)
               q <- quantile(s, probs = p, names = FALSE, type = 1)
               x[s > q,]
           },
           "sum.given.sum.exceeds" = {
               stopifnot(length(p) == 1)
               s <- rowSums(x)
               q <- quantile(s, probs = p, names = FALSE, type = 1)
               s[s > q,]
           }, stop("Wrong 'method'"))
}
