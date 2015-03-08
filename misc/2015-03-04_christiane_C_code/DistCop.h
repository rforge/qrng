#ifndef _Dist_H
#define _Dist_H

#define P0 (-3.22232431088E-1)
#define Q0 9.9348462606E-2
#define P1 (-1.0)
#define Q1 5.88581570495E-1
#define P2 (-3.42242088547E-1)
#define Q2 5.31103462366E-1
#define P3 (-2.04231210245E-2)
#define Q3 1.0353775285E-1
#define P4 (-4.53642210148E-5)
#define Q4 3.8560700634E-3

#define a0 2.50662823884
#define a1 -18.61500062529
#define a2 41.39119773534
#define a3 -25.44106049637
#define b0 -8.47351093090
#define b1 23.08336743743
#define b2 -21.06224101826
#define b3 3.13082909833
#define c0 0.3374754822726147
#define c1 0.9761690190917186
#define c2 0.1607979714918209
#define c3 2.76438810333863E-2
#define c4 3.8405729373609E-3
#define c5 3.951896511919E-4
#define c6 3.21767881768E-5
#define c7 2.888167364E-7
#define c8 3.960315187E-7

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

double MRG32k3a();

void InitMRG32k3a();

double Expon(double mu, double unif01);

double FDistNormale(double bidon, double x);

double InvNormalDist (double U);

double Normal(double Mean, double StdDeviation, double unif);

float VChisq(float df, float unif);

double *MOClayton(double theta, double *point, int dim);

double *CDMClayton(double theta, double *point, int dim);
#endif
