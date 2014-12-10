#include "rankone.h"
#include <stdlib.h>
#include <malloc.h>

static double **PointSet = NULL;

double **GenRank1 (int s, int n, int *GeneratingVector){
  int i,j;
  double *GenVectorOverN;

  GenVectorOverN = (double *) malloc (sizeof(double) * s);
  
  PointSet = malloc (n*sizeof(double*));

  for(i=0;i<n;i++){
    PointSet[i] = malloc (sizeof(double) *s);
  }

   
  
  for(j=0;j<s;j++){
    GenVectorOverN[j] = GeneratingVector[j]/((double) n);
    PointSet[0][j]=0;
  }
  printf("%d \n",GeneratingVector[0]);
  
  //start at i=1 since first point should be 0

  for (i=1;i<n;i++){
    for (j=0;j<s;j++){
      PointSet[i][j] = PointSet[i-1][j]+GenVectorOverN[j];
      if (PointSet[i][j] > 1) PointSet[i][j]=PointSet[i][j]-1.0;
    }
  }
  return(PointSet);
}

double **GenKorobovPrimElem(int s, int n, int a){
  
  int *GeneratingVector;
  int lastTerm;
  int i,j;
 
  PointSet = malloc (sizeof(double*)*n);
  GeneratingVector = malloc (sizeof(int)*s);
  for(i=0;i<n;i++){
    PointSet[i] = malloc (sizeof(double) *s);
  }
  GeneratingVector[0]=1;
  PointSet[1][0]=1.0/n;
  PointSet[0][0]=0;
  for(j=1;j<s;j++){
    GeneratingVector[j]=(GeneratingVector[j-1]*a)%n; 
    PointSet[1][j] = GeneratingVector[j]/((double) n);
    PointSet[0][j]=0;
  }
  lastTerm = GeneratingVector[s-1];
  for(i=2;i<n;i++){
    for(j=0;j<s-1;j++){
      PointSet[i][j] = PointSet[i-1][j+1];
    }
    lastTerm = (lastTerm * a) %n;
    if ((lastTerm==1) && (i<n-s+1)) printf("generator is not a primitive element\n");
    PointSet[i][s-1] = lastTerm/((double) n);
  }
  return(PointSet);
}
