#ifndef FOO_H
#define FOO_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP myCfun(SEXP arg1, SEXP arg2);

/**
 * C API
 */
double myCfun(double *arg1, const double arg2[]);

#endif