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

double Expon(double mu, double unif01);

double FDistNormale(double bidon, double x);

double InvNormalDist (double U);

double Normal(double Mean, double StdDeviation, double unif);

#endif
