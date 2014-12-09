/* FILE: qmc.h

   Header file for the QMC library.

   Mik Cieslak and Kris Luttmer | July 03, 2001 */


#ifndef _QMC_H
#define _QMC_H


/* includes for the QMC library */
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>


/*  defines for the QMC methods */
#define KOROBOV_METHOD    2
#define HALTON_METHOD     10
#define RHALTON_METHOD     18



/* defines for the InitQMC parameter list */

#define QMC_STATIC_PARAMS	7

#define QMC_METHOD 	0
#define QMC_DIMEN  	1
#define QMC_NPOINTS 2
#define QMC_RMETHOD 3
#define QMC_RANDS 4
#define QMC_RBASE   5
#define QMC_NPARAMS	6

/* defines for randomization methods */
#define NOSHIFT 1
#define ADDSHIFT 2
#define XORSHIFT 3			/* also known as the digital shift */

typedef unsigned int UINT; /* unsigned integer type define */

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------- */

int InitQMC ();
/* initializes the correct method, parameters described in qmc.c */

double *(*QMC) (void);
/* returns a point from the point set */

void (*ResetQMC) (void);
/* resets the qmc method */



void (*FreeMethod) (void);
/* frees the qmc method */

void FreeQMC();
/* -------------------------------------------------------------------- */

int LoadInputQMC(char * fname);

double GetQMC(int request);

int SetQMC(int request, double value);
/* -------------------------------------------------------------------- */

#define RANDVARONE 0
#define RANDVARTWO 1




/* -------------------------------------------------------------------- */

/* -------------------------------------------------------------------- */

double *GenRandom (void);
/* generate a random number to be used by AppRandom */

double *(*AppRandom) (double *point, double *r);
/* apply the randomization to point, and return the new point */

double RPoint (int currdim);
/* returns the a value for the (currdim) dimension from one point
   according to lastpoint and rpoint */



/* -------------------------------------------------------------------- */


#ifdef __cplusplus
};
#endif


#endif
