#ifndef _FUNCTIONQMC_H_
#define _FUNCTIONQMC_H_

//-----------------------------------------------------------------

double Hellekalek(double alpha,float dim,double *point);
double Sobol1994(float dim,double *point);
double Sobol2001( double *A,float dim,double *point);
double SobolAsot(double *A, float dim, double *point);
double SobolOld(float dim, double *point);
double Sobol2003(double c, float dim, double *point);
double OwenQMC(float dim,double *point);
double ARnold1(float dim,double *point);
double ARnold2(float dim,double *point);
double ARnold3(float dim,double *point);
double GenzQMC_PPeak(double *A,double *U,float dim,double *point);
double GenzQMC_Gaussian(double *A,double *U,float dim,double *point);
double GenzQMC_10(double *A,double *U,float dim,double *point);
double GenzQMC_Dis(double *A,double *U,float dim,double *point);
double KeisterQMC(float dim,double *point);
double PalphaPLR(float dim, double *point);
double PalphaPLRPlus1(int index, double *point);
double Bebe(int dim, double *point);
//------------------------------------------------------------------

#endif



