#include "rankone.h"
#include <stdio.h>
#include <malloc.h>


// ----------------------------------------------------------------------

int main (int argc, char **argv)
{  
  int i, j;
  int s,n,a;
  
 
 
  double **PointSet;
  int *Genvector;

  n = 11;
  a = 6;
  s = 10;
  
  Genvector = malloc (sizeof(int)*s);
  
  PointSet = malloc (n*sizeof(double*));

  for(i=0;i<n;i++){
    PointSet[i] = malloc (sizeof(double) *s);
  }

  Genvector[0]=1;

  for(j=1;j<s;j++){
    Genvector[j] = (Genvector[j-1]*a) % n;
  }
  PointSet = GenKorobovPrimElem(s,n,a);
  
  for(i=0;i<n;i++){
    for(j=0;j<s;j++){
      printf("%f ",PointSet[i][j]);
    }
    printf("\n");
  }
  PointSet = GenRank1(s,n,Genvector);
  for(i=0;i<n;i++){
    for(j=0;j<s;j++){
      printf("%f ",PointSet[i][j]);
    }
    printf("\n");
  }


  
  return (1);
}

// ----------------------------------------------------------------------
