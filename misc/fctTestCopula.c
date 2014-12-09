#include "qmc.h"
#include "FunctionQMC.h"
#include "DistCop.h"
#include <stdio.h>
#include <time.h>

// done in maple, see /home/clemieux/prog/testfct/keister.gnumeric


#define keis10 -0.5038691734
#define keis20 -0.8837024694
#define keis30 -0.6776504849
#define keis40 -0.2359661149
#define keis50 0.2275118812
#define keis120 0.1104770039
#define keis150 -0.6274167693
#define keis360 0.5891202116

//#define keisint 0.1104770039


#define NUMBLOCKS 1100
#define INTLOOP 0
#define EXTLOOP 1
#define INTERVAL 1000
#define c 1.0


// ----------------------------------------------------------------------

int main (int argc, char **argv)
{  
  int i, j, l;
  float bond;
  int k, m, n, nbl;
  
  double *point = NULL;
  double *A = NULL;
  double *rvector = NULL;
  double *rpoint = NULL;
  double result = 0.0;
  double var,mu,mu2,aver,prod;
  int s, realdim;
  double *cpoint;
  double *estim;

  clock_t start, end;
  double cpu_time_used;
  int CDM = 0;
  double theta = 1.0;

  
  

  if(!LoadInputQMC(argv[1]))
    {
      printf ("Cannot load the input file.\n");
      return (0);
    }
  
  
  
  if (!InitQMC ())
    {
      printf ("Cannot initialize the Quasi Monte Carlo method.\n");
      return (0);
    }
  
 
                                  
 
  
  s = (int) GetQMC(QMC_DIMEN);
  m = (int) GetQMC(QMC_RANDS);
  n = (int) GetQMC(QMC_NPOINTS);
  
  if (CDM==1) realdim = s; 
  else realdim = s-1;

  A = (double *) malloc (sizeof(double) * realdim);
                                
  cpoint = (double *) malloc (sizeof(double) * realdim);

  estim = (double *) malloc (sizeof(double) * m);

  for(i = 0; i< realdim;i++)
    {
      //A[i] = (nbvar-i)*(nbvar-i);
      A[i] = i+1;
      //A[i]=0.0;
      //A[i] = 0.01;
      //A[i] = (i+1)*(i+1);
    }  


  start = clock();

  for (i = 0; i < m; i++)
    {
      
      ResetQMC();

      rvector = GenRandom();
      aver = 0;      
      for (j = 0; j < n; j++)
	{
	  if ((point = QMC()) == NULL)
	    {
	      printf ("ERROR cannot generate point %d.\n", j);
	      FreeQMC();
	     
	      
	      return (0);
	    }
  
	  
	  if ((rpoint = AppRandom (point, rvector)) == NULL)
	    {
	      printf ("Cannot apply the randomization.\n");
	      FreeQMC();
	     
	      
	      return (0);
	    }


          cpoint = MOClayton(theta,rpoint,realdim);
          result = Sobol2003 (c, realdim, cpoint); 
            //result += KeisterQMC(nbvar,prpoint)/keisint;
	  //result = SobolAsot (A,realdim, cpoint); 
          aver += result;
	 
         
    
	}
      
      estim[i] = aver/n;


    }
  end = clock();
  mu=0;
  mu2=0;

  for(i=0;i<m;i++){
    mu += estim[i]/m;
    mu2 += estim[i]*estim[i]/(m-1);
  }

  var=mu2-mu*mu*m/(m-1);

  printf("rmethod = %d\n", (int) GetQMC(QMC_RMETHOD));
  printf("dimension = %d\n", s);
 
  printf("nb rand=%d\n", m);
  printf("n = %d\n", n);                                                                                
  printf("average value = %f\n",mu);
  printf("estimated variance = %f\n", var);    
  
  printf ("CPU time used: %f\n", cpu_time_used);

  FreeQMC();
 
  
  return (1);
}

// ----------------------------------------------------------------------
