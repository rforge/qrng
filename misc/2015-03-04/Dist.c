#include <math.h>
#include "Dist.h"

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

#define LIMITE 1E-16

double Expon(double mu, double unif01){
return -mu * log(1.0 - unif01);
} 


/*----------------------------------------------------------------*/
// from pages 95-96 of Kennedy and Gentle, "Statistical Computing"
// Dekker, New York, 1980.

double FDistNormale(double bidon, double x){

   double d1,d2,d3,d4,d5,d6,rep;
   int neg ;

   neg = 0;
   if (x<-64.0) return 0.0;
   else 
      if (x>64.0) return 1.0;
   else{
      d1= 0.0498673470;
      d2=0.0211410061;
      d3=0.0032776263;
      d4=0.0000380036;
      d5=0.0000488906;
      d6=0.0000053830;

      if (x<0.0){
         x = -x;
         neg = 1; 
      }
      rep=(1.0+d1*x+d2*x*x+d3*x*x*x+d4*x*x*x*x+d5*x*x*x*x*x+
			d6*x*x*x*x*x*x);
      rep=1.0-0.5*pow(rep,-16.0);

      if (neg==1) return 1.0-rep; else return rep;
   }
}

double InvNormalDist (double U)
{
  double Z;
  double Y;
  if (U>0.5) Y = sqrt(-(2.0*log(1.0-U)));
  else Y = sqrt(-(2.0*log(U)));
  Z = Y+((((Y*P4+P3)*Y+P2)*Y+P1)*Y+P0)
       /((((Y*Q4+Q3)*Y+Q2)*Y+Q1)*Y+Q0);
  if (U<0.5) Z = -Z;
  return Z;
} /* end InvNormalDist() */

/*------------------------------------------------------------------------*/
/* returns a normal r.v. of mean Mean, and standard deviation
StdDeviation, by transforming by inversion the uniform
0,1 unif*/

/*
double Normal(double Mean, double StdDeviation, double unif){

   if (StdDeviation < 0.0){  
       printf("erreur\n");  
       return 0.0;
   } 
   else if (unif < LIMITE){
     printf("Normale a la limite inf.\n");
     return -8.22;
   }
   else if (unif > 1-LIMITE){
     printf("Normale a la limite sup.\n");
     return 8.22;
   }
   else
      return Mean + InvNormalDist (unif) * StdDeviation;
}
*/

// Moros' algorithm, as given in Glasserman's book

double Normal(double Mean, double StdDeviation, double unif){
  double y, ay, r, x;
  if (StdDeviation < 0.0 || unif < 0.0 || unif > 1.0){  
       printf("erreur\n");  
       return 0.0;
   } 
  y = unif - 0.5;
  ay = fabs(y);
  if (ay < 0.42){
    r = y*y;
    x = y * (((a3 * r + a2) * r + a1 )*r + a0)/((((b3 * r + b2)*r + b1)*r + b0)*r +1);
  }
  else{
    if(unif==0) unif = 0.00000000045;
    r = unif;
    if (y > 0) r = 1- unif;
    r = log(-logl(r));
    x = c0 + r*(c1+r*(c2+r*(c3+r*(c4+r * (c5+r*(c6+r*(c7+r*c8)))))));
    if(y<0) x = -x;
  }
  //if(fabs(x)>20) printf("unif=%15.10f\n",unif);
  return(Mean + StdDeviation* x);
}





