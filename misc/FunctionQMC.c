#include <math.h>
#include <stdio.h>
#include "DistCop.h"
#include "FunctionQMC.h"


#define PIE 3.14159265358979


//----------------------------------------------------------------------

double Hellekalek(double alpha, float dim, double *point)
{
  int i;
  double temp1, temp2;
  double result;
  
  result = 1;
  temp2 = 1/(alpha + 1);

  
  for (i = 0; i < dim; i++){
    
    temp1 = pow(point[i], alpha);
    
    result = result * (temp1 - temp2);
    
  }
      
      
  
    return(result);


}//end of Hellekalek
//--------------------------------------------------------------------------

double Sobol1994(float dim, double *point)
     
{
  int i;
  double temp1, temp2, result;
  temp1 = 0;
  temp2 = 0;
  result = 1;
 
  for(i = 0; i < dim; i++){
     
    temp1 = i + 2 * point[i];
    temp2 = i + 1;
    result = result * (temp1 / temp2);

   
  }


  return(result);
  }//end of sobolQMC
//-------------------------------------------------------------------------

double Sobol2001(double *A, float dim, double *point)
     
{
  int i;
  double temp1, temp2, result, value;
  result = 1; 
 
  for(i = 0; i < dim; i++){


    if(i <= dim/4){
      A[i] = 0;
      temp1 = fabs(4 * point[i] - 2)+ A[i];
      temp2 = 1 + A[i];
    }


    if (i > dim/4){
      A[i] = 3;
      temp1 = fabs(4 * point[i] - 2)+ A[i];
      temp2 = 1 + A[i];
    }


    
    result = result * (temp1/temp2);
  }
  
    return(result);
}//end of Sobol

//-------------------------------------------------------------------------

double SobolAsot(double *A, float dim, double *point)

{
  int i;
  double temp1, temp2, result, value;
  result = 1;

  for(i = 0; i < dim; i++){
 
      temp1 = fabs(4 * point[i] - 2)+ A[i];
      temp2 = 1 + A[i];
      result = result * temp1/temp2;
    }
  return (result);
}
//-------------------------------------------------------------------------

double SobolOld(float dim, double *point)
     
{
  int i;
  double temp1, temp2, result, value;
  result = 1.0; 
 
  for(i = 0; i < dim; i++){


      result *= fabs(4 * point[i] - 2.0);
   
   
  }
  
    return(result);
}//end of Sobol

//-------------------------------------------------------------------------

double Sobol2003(double c, float dim, double *point)
     
{
  int i;
  double result;
  result = 1.0; 
 
  for(i = 0; i < dim; i++){


      result *= 1.0 + c *(point[i] - 0.5);
   
  }
  
 return(result);
}//end of Sobol
//---------------------------------------------------------------------------

double OwenQMC(float dim, double *point)
{
  int i;
  double temp1, temp2, result;
  result = 1;

  for(i = 0; i < dim; i++){
    temp1 = point[i] - 0.5;
    result = result * temp1;
  }


  temp2 = pow(12, (dim/2.0));
  result = result * temp2;

  return(result);
}//end of OwenQMC
//---------------------------------------------------------------------------

double ARnold1(float dim, double *point)
{
  int i;
  double temp1, temp2, result;
  temp1 = 0;


  for (i = 0; i < dim; i++)
    {
      temp1 = temp1 + fabs(4 * point[i] - 2);
    }


  temp2 = 1/dim;
  result = temp2 * temp1;
  return(result);

}//end Arnold1y
//--------------------------------------------------------------------------

double ARnold2(float dim, double *point)
{
  int i;
  double result;
  result = 1;
 

  for(i = 0; i < dim; i++)
    {
      result =  result * fabs(4 * point[i] - 2); 



    }
  return(result);
}//end ARnold2
//--------------------------------------------------------------------------

double ARnold3(float dim, double *point)
{
  int i;
  double temp, result;
  result = 1;


  for(i = 0; i < dim; i++)
    {
      result = result * ((PIE/2) * sin(PIE * point[i]));
    }


  return(result);
}//end ARNOLD3

//-------------------------------------------------------------------------

double GenzQMC_PPeak(double *A, double *U, float dim, double *point)
{
  int i;
  double temp1, result;
  result = 1;


  for(i = 0; i < dim; i++)
    {
      temp1 = 1/( pow(A[i], -2) + pow((point[i] - U[i]), 2));
      result = result * temp1;
    }


  return(result);
}//end of GenzQMC1
//---------------------------------------------------------------------------

double GenzQMC_Gaussian(double *A, double *U, float dim, double *point)
{
  int i;
  double temp1, result; 
  temp1 = 0;
  

 for (i = 0; i < dim; i++)
   {
      temp1 = temp1 + (pow(A[i], 2) * pow((point[i] - U[i]), 2));
   }


 result = exp( -temp1);

 
 return(result);
}//end of  GenzQMC_Gaussian
//--------------------------------------------------------------------------

double GenzQMC_10(double *A,  double *U,  float dim,  double *point) 
{
  int i;
  double temp1, result;
  temp1 = 0;

  for(i = 0; i < dim; i++)
    {
      temp1 = temp1 + ( A[i] * fabs(point[i] - U[i]));
    }


  result = exp(-temp1);
  
  return(result);
}//end of GenzQMC_10
//-------------------------------------------------------------------------

double GenzQMC_Dis(double *A,double *U,float dim,double *point)
{
  int i;
  double temp,result;
  temp = 0;
  

  if ((point[1] < U[1])  ||  (point[2] < U[2])) 
      return(0);
  else
    {

      for(i = 0; i < dim; i++)
	{
	  temp = temp +( A[i] * point[i]);
	}


      result = exp(-temp);
      return(result);
    }
}//end of GenzQMC_Dis
//--------------------------------------------------------------------------
double KeisterQMC(float dim,double *point)
{
  // we don't multiply by pi^(s/2), compared to Keister and, e.g., Hong and Hickernell

  int i;
  double temp,result,u,con,fac;
  temp = 0;
 
  con = sqrt(PIE);
  fac = 1;

  for(i = 0; i < dim; i++) 
    {
      u = Normal (0.0, 1.0, point[i]);
      temp = temp + u*u ; 
      //fac = fac* con;
    }

  temp = temp * 0.5;
  //if (temp<-20 || temp > 20) printf("temp=%f\n",temp);
  temp = cos(sqrt(temp));
  //printf("%f\n",temp);
  
  result = temp;

  return (result);
}//end of KeisterQMC

double PalphaPLR(float dim, double *point){

  int i;
  double value, temp;

  temp = 1;
  for(i = 0; i < dim; i++){
    value = log(point[i])/ log(2);
    value = floor(value);
    temp = temp * (3 - 6 * pow(2.0, value));
    
  }
  temp = temp -1;
  return (temp);
}

double PalphaPLRPlus1(int index, double *point){

  int i;
  double value;

  //printf("%f\n", point[index]);
  value = log(point[index])/ log(2);
  value = floor(value);
 
  return (3 - 6 * pow(2.0, value));
}
//--------------------------------------------------------------------------



double Bebe(int dim, double *point){
  int j;
  double sum = 0.0;

  for(j=0;j<dim;j++){
    sum += point[j];
  }
  return(sum);
}





