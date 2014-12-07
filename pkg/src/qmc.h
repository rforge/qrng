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
#define MONTECARLO_METHOD 1
#define KOROBOV_METHOD    2
#define PKOROBOV_METHOD	  3	
#define SOBOL_METHOD      4
#define SHIFTNET_METHOD   5
#define GFAURETT94_METHOD 6
#define GFAUREF01_METHOD  7
#define FAURE_METHOD      8
#define GENGFAURE_METHOD      9
#define HALTON_METHOD     10
#define GHALTON_METHOD    11
#define SALZBURG_METHOD   12
#define LHS_METHOD        13
#define MLHS_METHOD       14
#define GENERIC_NET_METHOD   15
#define ILHS_METHOD   16
#define DIAFAURE_METHOD      17
#define RHALTON_METHOD     18
#define AHALTON_METHOD     19           // Atanassov + random shift
#define DAHALTON_METHOD     20          // Atanassov with no random shift
#define KWHALTON_METHOD     21  // deterministic Kocis-Whiten
#define RKWHALTON_METHOD     22  // random Kocis-Whiten
#define DIAFAURERAND_METHOD      23 // DET DIa + random lower part of NLT
#define PERMHALTON_METHOD      24 // random permutations like Okten
#define RPERMHALTON_METHOD      25 // random permutations like Okten
#define SCRAHALTON_METHOD      26 // AD but with det. scrambling based on power of admissible factors; same element on diagonal and subdiagonal; also random shift
#define BADSOBOL_METHOD      27 // 


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
#define OWENSCRAMBLE 4
#define PNLTSCRAMBLE 5 /* NLT Scramble with coefficients in BF(k) for some k = 10th param */
#define PDIASCRAMBLE 6 /* DIA Scramble with coefficients in BF(k) for some k = 10th param */
#define PNLTSCRAMBLEEXCDIA 7

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

void (*ResetQMC2) (void);
/* resets the qmc method but just the points (not the matrices in Faure)*/


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

typedef struct
{
  double average[2]; 	/* average of var 1 */
  double squares[2]; 	/* average of squares of var 1 */
  double totalaverage; 	/* average of both var 1 and 2 */
  double observations; 	/* num of observations */
} BLOCKQMC;

/* define a STATSQMC type */
typedef struct
{
  int numblocks;
  BLOCKQMC *blocks;
} STATSTYPE;

typedef STATSTYPE* STATSQMC;

/* define a CONFIDENCE INTERVAL type */
typedef struct
{
  double low;
  double high;
} CONFIDENCEINTERVAL;


/* -------------------------------------------------------------------- */

STATSQMC InitStat (int numblocks);
/* initializes and returns a STATSQMC type with blocks, numblocks */

void StatUpdate1 (STATSQMC stats, int block, double entry);
void StatUpdate2 (STATSQMC stats, int block, double entry1, double entry2);
/* updates randvar with entry at block with index block */

double StatAverage (STATSQMC stats, int block, int randvar);
/* returns average for variable randvar */

double StatVariance (STATSQMC stats, int block, int randvar);
/* returns variance for variable randvar */

CONFIDENCEINTERVAL StatConfIntval (STATSQMC stats, int block, int randvar,
                                   double level);
/* returns confidence interval for variable randvar */

double StatCovariance (STATSQMC stats, int block);
/* returns the covariance of the two random variables at block, block */

double StatCorrelation (STATSQMC stats, int block);
/* returns the correlation of the two random variables at block, block */

void ResetStat (STATSQMC stats, int block);
/* reset the stats for block at index block */

void FreeStat (STATSQMC stats);
/* frees all of the blocks in stats */

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
