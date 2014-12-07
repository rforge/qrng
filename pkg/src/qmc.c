// FILE: qmc.c

/* Implementation of the QMC Library. 
   See the Header file qmc.h for a full description of the functions.
   Only one function is implemented, InitQMC, because the rest of the
   functions are just function pointers.  InitQMC assigns those
   functions the appropriate QMC method funtions.
   The parameters for InitQMC is an unsigned integer array of size
   MAXARGS.  There are defines in qmc.h that describe the contents
   of each index. 

 Kris Luttmer and Mik Cieslak | March 11, 2001
 Updated February 25th,  2002
*/

#include "qmc.h"
#include "lhs.h"
#include "sal.h"
#include "halton.h"
#include "rhalton.h"
#include "ahalton.h"
#include "scrahalton.h"
#include "dahalton.h"
//#include "kwhalton.h"
#include "permhalton.h"
#include "rpermhalton.h"
//#include "rkwhalton.h"
#include "shift.h"
#include "mc.h"
#include "faure.h"
#include "korobov.h"
#include "sobol.h"
#include "badsobol.h"
#include "combinedpk.h"
#include "ghalton.h"
#include "randqmc.h"


static double * args = NULL;


// -----------------------------------------------------------------------
int InitQMC ()
/* Purpose: to prepare the correct QMC method for generating points.
       Pre: parameters[MAXARGS] must contain correct values as described
            by the defines for it in qmc.h.
      Post: the choosen method is initialized with the correct parameters.
            the three function pointers are set to the correct methods. 
            Upon success 1 is returned else 0. */
{
  int retvalue; /* the return value */
 
  switch ((int)args[QMC_METHOD])
    {
    case MONTECARLO_METHOD:
     
      retvalue = InitMC(args[QMC_DIMEN]);
      ResetQMC = ResetMC;
      FreeMethod = FreeMC;
      QMC = MC;      
		break;
    
    case KOROBOV_METHOD:
      retvalue = InitKorobov (args[QMC_DIMEN], args[QMC_NPOINTS],
                             args[QMC_NPARAMS+1], args[QMC_NPARAMS+2]);
      QMC = Korobov;
      ResetQMC = ResetKorobov;
      FreeMethod = FreeKorobov;
      break;
    
    case PKOROBOV_METHOD:
      retvalue = InitCpKorobov (args[QMC_DIMEN], args[QMC_NPOINTS],
                             args[QMC_NPARAMS+1], &args[QMC_NPARAMS+2]);
      QMC = CpKorobov;
      ResetQMC = ResetCpKorobov;
      FreeMethod = FreeCpKorobov;
      break;
  
    case SOBOL_METHOD:
      retvalue = InitSobol (args[QMC_DIMEN], args[QMC_NPOINTS], args[QMC_NPARAMS+1]);
      QMC = Sobol;
      ResetQMC = ResetSobol;
      FreeMethod = FreeSobol;
      break;

    case BADSOBOL_METHOD:
      retvalue = InitBadSobol (args[QMC_DIMEN], args[QMC_NPOINTS], args[QMC_NPARAMS+1]);
      QMC = BadSobol;
      ResetQMC = ResetBadSobol;
      FreeMethod = FreeBadSobol;
      break;
    
    case SHIFTNET_METHOD:
      retvalue = InitShiftNet (args[QMC_DIMEN], args[QMC_NPOINTS]);
      QMC = ShiftNet;
      FreeMethod = FreeShiftNet;
      ResetQMC = ResetShiftNet;
      break;
   
    case GFAURETT94_METHOD:
    case GFAUREF01_METHOD:
    case FAURE_METHOD:
    case DIAFAURE_METHOD:
    case DIAFAURERAND_METHOD:
    case GENGFAURE_METHOD:
      //printf("method is %f\n",args[QMC_METHOD]);
      retvalue = InitFaure (args[QMC_DIMEN], args[QMC_NPARAMS+1],
                            args[QMC_METHOD], args[QMC_NPOINTS], 
                            args[QMC_NPARAMS+2], args[QMC_NPARAMS+3]);
      QMC = Faure;
      FreeMethod = FreeFaure; 
      ResetQMC = ResetFaure;
      ResetQMC2 = ResetFaurePts;
      break;
    case HALTON_METHOD:
      retvalue = InitHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1], args[QMC_NPARAMS+2]);
      if (args[QMC_NPARAMS+1]<8)
        QMC = Halton;
      else{
	QMC = FaurePermHalton;
        printf("got it");
      }
      FreeMethod = FreeHalton;
      ResetQMC = ResetHalton;
      break;


    case GHALTON_METHOD:
      retvalue = InitGHalton (args[QMC_DIMEN]);
      QMC = GHalton;
      FreeMethod = FreeGHalton;
      ResetQMC = ResetGHalton;
      break;

    case RHALTON_METHOD:
      retvalue = InitRHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1], args[QMC_NPARAMS+2]);
      if( args[QMC_NPARAMS+1] <8)
        QMC = RHalton;
      else
        QMC = RFaurePermHalton;
      FreeMethod = FreeRHalton;
      ResetQMC = ResetRHalton;
      break;

    case AHALTON_METHOD:
      retvalue = InitAHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = AHalton;
      FreeMethod = FreeAHalton;
      ResetQMC = ResetAHalton;
      break;

    case SCRAHALTON_METHOD:
      retvalue = InitScrAHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = ScrAHalton;
      FreeMethod = FreeScrAHalton;
      ResetQMC = ResetScrAHalton;
      break;

    case DAHALTON_METHOD:
      retvalue = InitDAHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = DAHalton;
      FreeMethod = FreeDAHalton;
      ResetQMC = ResetDAHalton;
      break;
      /*
    case KWHALTON_METHOD:
      retvalue = InitKWHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = KWHalton;
      FreeMethod = FreeKWHalton;
      ResetQMC = ResetKWHalton;
      break;

    case RKWHALTON_METHOD:
      retvalue = InitRKWHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = RKWHalton;
      FreeMethod = FreeRKWHalton;
      ResetQMC = ResetRKWHalton;
      break;*/ 
    case PERMHALTON_METHOD:
      retvalue = InitPermHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = PermHalton;
      FreeMethod = FreePermHalton;
      ResetQMC = ResetPermHalton;
      break;
    case RPERMHALTON_METHOD:
      retvalue = InitRPermHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1]);
      QMC = RPermHalton;
      FreeMethod = FreeRPermHalton;
      ResetQMC = ResetRPermHalton;
      break;
    case SALZBURG_METHOD:
      retvalue = InitSalzburg(args[QMC_DIMEN],args[QMC_NPARAMS+2], args[QMC_NPARAMS+1], 
	                          args[QMC_NPARAMS+3], args[QMC_NPARAMS+4]);
      QMC = Salzburg;
      FreeMethod = FreeSalzburg;
      ResetQMC = ResetSalzburg;
      break;

    case GENERIC_NET_METHOD:
      retvalue = InitGenericNet(args[QMC_DIMEN],args[QMC_NPARAMS+2], args[QMC_NPARAMS+1]);
      if (args[QMC_NPARAMS+2] == 2) 
        QMC = SalzburgBase2;
      else
        QMC = Salzburg;
      FreeMethod = FreeSalzburg;
      ResetQMC = ResetSalzburg;
      break;
        
    case LHS_METHOD:
      retvalue = InitLHS(args[QMC_DIMEN], args[QMC_NPOINTS]);
      QMC = LHS;
      FreeMethod = FreeLHS;
      ResetQMC = ResetLHS;
      break;
    
    case MLHS_METHOD:
      retvalue = InitLHS(args[QMC_DIMEN], args[QMC_NPOINTS]);
      QMC = MLHS;
      FreeMethod = FreeLHS;
      ResetQMC = ResetLHS;
      break;

    case ILHS_METHOD:
      retvalue = InitLHS(args[QMC_DIMEN], args[QMC_NPOINTS]);
      QMC = ILHS;
      FreeMethod = FreeLHS;
      ResetQMC = ResetLHS;
      break;

	 default:
      retvalue = 0;
      QMC = 0;
      ResetQMC = 0;
      FreeMethod = 0;
      break;
  }

  if(retvalue == 0)
 	return 0;
 
  retvalue = InitRandom(args);
  
  if(retvalue == 0)
  {
  	fprintf(stderr, "Cannot Initialize Random Method\n");
  	FreeMethod();
 	retvalue = 0;
  }
 
  return (retvalue);
}

// -----------------------------------------------------------------------
void FreeQMC()
{
	if(FreeMethod)
		FreeMethod();
		
	FreeRandom();
	
	if(args)
		free(args);
}

// -----------------------------------------------------------------------

int LoadInputQMC(char * fname)
{
	FILE * infile;
	double tempargs[QMC_STATIC_PARAMS];
	int i;
	
	infile = fopen(fname, "r");
	
	if(infile == NULL)
	{
		fprintf(stderr, "Cannot open input file %s\n", fname);
		return 0;
	}
	
	for(i = 0; i < QMC_STATIC_PARAMS; i++)
		fscanf(infile, "%lf", &tempargs[i]);
	
	args = malloc(sizeof(double)*(QMC_STATIC_PARAMS + tempargs[QMC_NPARAMS]));
	
	for(i = 0; i < QMC_STATIC_PARAMS; i++)
		args[i] = tempargs[i];
	
	for(i = 0; i < args[QMC_NPARAMS]; i++)
		fscanf(infile, "%lf", &args[i+QMC_STATIC_PARAMS]);
		
	return 1;
}
	
double GetQMC(int request)
{

	if(args == NULL)
	{
		fprintf(stderr, "Cannot get QMC parameter without loading input file first\n");
		return 0;
	}
	
	if(request < QMC_STATIC_PARAMS + args[QMC_NPARAMS])
		return args[request];
	
	fprintf(stderr, "QMC Parameter requested is out of range\n");
	return 0;
}

int SetQMC(int request, double value)
{
	if(args == NULL)
	{
		if(request != QMC_NPARAMS)
		{
			fprintf(stderr, "Must set the number of extra parameters first\n");
			return 0;
		}
		
		args = malloc(sizeof(double)*(value + QMC_STATIC_PARAMS));
		args[QMC_NPARAMS] = value;
		return 1;
	}
		
	
	if(request < QMC_STATIC_PARAMS + args[QMC_NPARAMS])
	{	
		args[request] = value;
		return 1;
	}	
	
	fprintf(stderr, "QMC Parameter requested is out of range\n");
	return 0;
}








