#include <malloc.h>
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

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

// these are the seeds: s1* must be < m1 and not all 0; 
//                      s2* must be < m2 and not all 0


double  s10, s11, s12, s20, s21, s22;

double MRG32k3a ()
    {
    long   k;
    double p1, p2;
    /* Component 1 */
    p1 = a12 * s11 - a13n * s10;
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10 = s11;   s11 = s12;   s12 = p1;
    /* Component 2 */
    p2 = a21 * s22 - a23n * s20;
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20 = s21;   s21 = s22;   s22 = p2;
    /* Combination */
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm);
    }

// end of PRNG

// ========================================================================
void InitMRG32k3a ()
{
  s10 = 12345, s11 = 56789, s12 = 321, s20 = 654, s21 = 987, s22 = 4321;
}



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
   else
      return Mean + InvNormalDist (unif) * StdDeviation;
}

*/
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
    r = unif;
    if (y > 0) r = 1- unif;
    r = log(-log(r));
    x = c0 + r*(c1+r*(c2+r*(c3+r*(c4+r * (c5+r*(c6+r*(c7+r*c8)))))));
    if(y<0) x = -x;
  }
  return(Mean + StdDeviation* x);
}


// From Bratley, Fox and Schrage "Guide to Simulation", 1987

float VChisq(float df, float unif){
	
   float sqp5 = 0.7071067811;
   float dwarf = 0.000000000000001;
   float z, vchisq, arg, sqdf, ch, zsq;

   if (df == 1.0){ 
      z = InvNormalDist(1.0-(1.0-unif)/2.0);
      vchisq = z * z;
      return vchisq;
   }
   else{
      if (df==2.0){
	 arg=1.0-unif;
	 if(arg<dwarf) arg=dwarf; 
	    vchisq=-2.0*log(arg);
	    return vchisq;
      }
      else{
         if(unif>dwarf){
	    z = InvNormalDist(unif);
	    sqdf = sqrt(df);
	    zsq = z * z;
	    ch=-(((3753.0*zsq+4353.0)*zsq-289517.0)*zsq-289717.0)*z*sqp5/9185400.0;
	    ch=ch/sqdf+(((12.0*zsq-243.0)*zsq-923.0)*zsq+1472.0)/25515.0;
	    ch=ch/sqdf+((9.0*zsq+256.0)*zsq-433.0)*z*sqp5/4860.0;
	    ch=ch/sqdf-((6.0*zsq+14.0)*zsq-32.0)/405.0;
	    ch=ch/sqdf+(zsq-7.0)*z*sqp5/9.0;
	    ch=ch/sqdf+2.0*(zsq-1.0)/3.0;
	    ch=ch/sqdf+z/sqp5;
	    vchisq=df*(ch/sqdf+1.0);
	    return vchisq;
	 }
         else return 0.0;
      }
   }
}

double *MOClayton(double theta, double *point, int dim)
{
  int i,status,mywhich;
  double bound;
  
  double unsurTheta,myp,myq;
  double *cpoint,*V;
  
  mywhich = 2;
 
  bound = -1.0;
  

  cpoint = (double *) malloc (sizeof(double) * dim);
  V = (double *) malloc (sizeof(double));

  
  unsurTheta = 1.0/theta;
  
  myp = point[dim];
  //////ORDERING///////////////////////
  //myp = point[0];
  
  myq = 1.0-myp; 
  gaminv(&unsurTheta,V,&bound,&myp,&myq,&status);
  //cdfgam(&mywhich,&myp,&myq, V, &unsurTheta,&myone,&status,&bound);
  //printf("V=%f\n",*V);
   
  //cpoint[0] = pow(1-log(point[0])/(*V),-unsurTheta);
  // phi = Normal (0.0, 1.0, cpoint[0]);
  
  for(i=0;i<dim;i++){
    
    //c = pow(point[i],-theta);
    //c2 = pow(point[i],-1.0/(i+unsurTheta)); // + 1.0 apres i?? 9 juillet
    cpoint[i] = pow(1-log(point[i])/(*V),-unsurTheta);
    //cpoint[i] = point[i]; // revert to indep
  }

  return(cpoint);
}

double *CDMClayton(double theta, double *point, int dim)
{
  int i;
  double c, sum2, S, phi,cc2;
  double treal, sum, deltat;
  double temp1, temp2,sqT,unsurTheta;
  double *cpoint;
  int *index;

  cpoint = (double *) malloc (sizeof(double) * dim);
  index = (int *) malloc (sizeof(int) * dim);

  for(i=0;i<dim;i++){
    //index[i]=dim-i-1;
    index[i]=i;
  }

  unsurTheta = 1.0/theta;
  //printf("unsurTheta=%f\n",unsurTheta);
  cpoint[0] = point[index[0]];
  c = pow(cpoint[0],-theta);
  //phi = Normal (0.0, 1.0, cpoint[0]);
  sum2 = 0; 
  for(i=1;i<dim;i++){
    sum2 += c;
    //printf("sum2=%f\n",sum2);
    
    cc2 = pow(point[index[i]],-1.0/(((double) i)+unsurTheta)); // + 1.0 apres i?? 9 juillet
    cpoint[i] = pow((1+(1-i+sum2)*(cc2-1.0)),-unsurTheta);
    c = pow(cpoint[i],-theta);
    //cpoint[i] = point[i]; // revert to indep
  }
  return(cpoint);
}



