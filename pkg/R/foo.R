myRfun <- function(arg1, arg2)
{
    stopifnot(arg1 > 0, arg2 > 0)
    .Call(myCfun, arg1, arg2)
}
