#include "qmc.h"
#include "AsianOption.h"
#include "Dist.h"

#define NUMBLOCKS 1100
#define INTLOOP 0
#define EXTLOOP 1
#define INTERVAL 1000

// ----------------------------------------------------------------------
int main (int argc, char **argv)
{  
  int i, j, m, n, s, l;
  int k;

  double sumalpha=0.0;
		
  double *point = NULL;
  double *rvector = NULL;
  double *rpoint = NULL;
  double *alpha = NULL;
  double result = 0.0;

  STATSQMC stats = NULL;
  CONFIDENCEINTERVAL cf;


  /*double K = 100;
  double T1 = 0.246576;
  double S0 = 100.0;
  double r = 0.086177;
  double T = 0.328767;
  double sigma = 0.2;*/

  double K = 50;
  //double realvalue = 2.105612;
  //double realvalue = 7.047734;
  //double realvalue= 4.052543; // s40
  //double realvalue= 4.015592;  // s75 K=50
  //double realvalue = 7.015250; // s75 K = 45
  //double realvalue = 2.073027; // s75 K = 55
  //double realvalue = 1.003532; // s =40, K = 60
  //double realvalue = 26.392012; // dig, S0 =50
    //double realvalue = 0;
  double T1 = 0.0;
  double S0 = 50.0;
  double r = 0.05;
  double T = 1.0;
  double theta = 6.0;
  //double T = 2.4615385;
  //double T = 0.0769;
  double sigma = 0.3;
  double salpha = 1.0;
  int nbl,realdim;
  int CDM = 0;

  InitMRG32k3a();

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

  m = GetQMC(QMC_RANDS);
  n = GetQMC(QMC_NPOINTS);
  s = GetQMC(QMC_DIMEN);

  if ((stats = InitStat (NUMBLOCKS)) == NULL)
  {
     printf ("Cannot initialize stats.\n");
     FreeQMC();
     return (0);
  }

  

  if (CDM==1) realdim=s; else realdim = s-1;
  alpha = (double *) malloc (sizeof(double) * realdim);
  
  printf("S0=%f,theta=%f,s=%d\n",S0,theta,s);

  if (salpha==1.0){
   for(j=0;j<realdim;j++){
     alpha[j] =1.0;
    }
  }

  else{

    sumalpha = (1.0-pow(salpha,realdim))/(1.0-salpha);

    alpha[0]=salpha;

    for(j=1;j<realdim;j++){
      alpha[j] = salpha*alpha[j-1];
    }

    for(j=0;j<realdim;j++){
      alpha[j] = realdim*alpha[j]/sumalpha;
    }
  }
  
  //printf("gecko are awesome\n");

  for (i = 0; i < m; i++)
  {
     ResetStat (stats, INTLOOP);
     ResetQMC();

     rvector = GenRandom();
     //for(j=0;j<s;j++) rvector[j]=0;
     //printf("rvector[0]=%f\n",rvector[0]); 
     nbl = 2;
     for (j = 0; j < n; j++)
     {
       if ((point = QMC()) == NULL)
       {
	  printf ("ERROR cannot generate point %d.\n", j);
	  FreeQMC();
	  FreeStat (stats);
	  return (0);
       }
       
       

       if ((rpoint = AppRandom (point, rvector)) == NULL)
       {
	  printf ("Cannot apply the randomization.\n");
	  FreeQMC();
	  FreeStat (stats);
	  return (0);
       }
       //if ((j<3)&&(j>0) && i < 5)
       //printf("%f %f %f\n",point[0],point[9], point[19]);
       //if (j<204) {
		//for(k=0;k<3;k++) printf("%15.10f, ",Normal(0,1,point[k]));
                //printf("%15.10f\n",Normal(0,1,point[3]));
		// }
       printf("%f %f\n",point[11],point[23]);
       //rpoint[realdim-1]=MRG32k3a();
       //result = BaskOptionClayton (K, S0, r, T, sigma, realdim, theta, alpha, rpoint);           
       //printf("%d,%f\n",s-1,rvector[s-1]);
       result = BaskOptionClaytonMO (K, S0, r, T, sigma, realdim, theta, alpha, rpoint);
       //result = BaskOptionIndep (K, S0, r, T, sigma, s, rpoint);
       //result = DigitalOption (S0, r, T1, T, sigma, s, rpoint);
       //printf("res=%f\n",result);
       StatUpdate1 (stats, INTLOOP, result);
      
       if (j >0 && (((j+1) % INTERVAL) == 0)){
	 StatUpdate1 (stats, nbl, StatAverage (stats, INTLOOP, RANDVARONE));
         //StatUpdate1 (stats, nbl+1, fabs(realvalue - StatAverage(stats,INTLOOP,RANDVARONE)));
         if (i == (m -1)){
	   printf("%15.10f\n", StatAverage(stats, nbl, RANDVARONE));
           printf("%15.10f\n",StatVariance(stats, nbl, RANDVARONE));
	 }
         nbl=nbl+2;
       }
       
     }
      
     result = StatAverage (stats, INTLOOP, RANDVARONE);
     StatUpdate1 (stats, EXTLOOP, result);
  }
 
  printf ("Average: %f\n", StatAverage (stats, EXTLOOP, RANDVARONE));
  printf ("Variance: %15.10e\n", StatVariance (stats, EXTLOOP, RANDVARONE));
  cf = StatConfIntval (stats, EXTLOOP, RANDVARONE, 0.95);
  printf ("Confidence Interval (95%%): %f to %f\n", cf.low, cf.high);
 

  FreeQMC();
  FreeStat (stats);
  return (1);
}

// ----------------------------------------------------------------------
