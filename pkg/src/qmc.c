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

#include "rhalton.h"

#include "korobov.h"

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
    
    case KOROBOV_METHOD:
      retvalue = InitKorobov (args[QMC_DIMEN], args[QMC_NPOINTS],
                             args[QMC_NPARAMS+1], args[QMC_NPARAMS+2]);
      QMC = Korobov;
      ResetQMC = ResetKorobov;
      FreeMethod = FreeKorobov;
      break;
    
    
    


    

    case RHALTON_METHOD:
      retvalue = InitRHalton (args[QMC_DIMEN],args[QMC_NPARAMS+1], args[QMC_NPARAMS+2]);
      
      QMC = RHalton; //args[QMC_NPARAMS+1] =1 or 6 depending on original or gene.
      
      FreeMethod = FreeRHalton;
      ResetQMC = ResetRHalton;
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








