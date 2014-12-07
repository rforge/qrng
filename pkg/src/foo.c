#include <Rmath.h>
#include "foo.h"

/**
 * Compute ...
 *
 * @param
 * @return
 * @author
 */
SEXP myCfun(SEXP arg1, SEXP arg2)
{
    int arg1 = PROTECT(arg1_ = asInteger(arg1));
    SEXP res = allocVector(REALSXP, arg1); // result
    // main part
    GetRNGstate();
    res = arg1 + arg2; // call Christiane's C functions
    PutRNGstate();
    // return
    UNPROTECT(1);
    return(res);
}