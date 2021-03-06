/* C function for computing a Sobol sequence **********************************/

#include "sobolC.h"


/**
 * @title Initialization of the Array v
 * @param v (sobolMaxDim, sobolMaxCol)-matrix
 * @param n number of points
 * @param d dimension
 * @return void
 * @author Marius Hofert based on C. Lemieux's RandQMC
 */
void init_array_v (unsigned int **v, int n, int d)
{
	int i, j, k;
	int newv, temp, degree, numcols;

	/* First the first dimension */
	numcols = ceil(log(n)/log(2));
	for (i = 0; i < numcols; i++)
		v[0][i] = 1;

	/* Now the rest of the dimensions */
	for (i = 1; i < d; i++)
	{
		/* Find the degree of polynomial i */
		degree = sobolMaxDegree;
		while (!((poly[i] >> degree) & 0x1))
			degree--;

		/* Copy known values from minit */
		memcpy(v[i], minit[i-1], degree * sizeof(unsigned int));

		/* Find the rest of the values */
		for (j = degree; j < numcols; j++)
		{
			newv = v[i][j-degree];
			temp = 1;
                        for (k = degree-1; k >= 0; k--) {
				if ((poly[i] >> k) & 0x1)
					newv = newv ^ (v[i][j-(degree-k)] << temp);
				temp++;
			}
			v[i][j] = newv;
		}
	}

	/* Calculate corresponding v<i,j>, v<i,j> = m<i,j>/2^j */
	temp = 1;
	for (j = numcols-2; j >= 0; j--){
		for (i = 0; i < d; i++)
			v[i][j] = v[i][j] << temp;
		temp++;
	}
}

/**
 * @title Generate n Points of a d-dimensional Sobol Sequence
 * @param n number of points
 * @param d dimension
 * @param randomize logical (here: integer) indicating whether a digital
 *        shift should be included
 * @param res pointer to the result matrix
 * @return void
 * @author Marius Hofert based on C. Lemieux's RandQMC
 */
void sobol(int n, int d, int randomize, double *res)
{
	int i, count, numcols;
	int column; /* the column to use */
	unsigned int maxn;
	double U; /* temporal storage for a single element of point */
	double recipd = 0.0;
	unsigned int *lastpoint = NULL;
	unsigned int randint, point; /* the int representation of U */
	int rmaxcoeff = (int) ((sizeof(int) * 8));
	double rmaxint = pow (2, rmaxcoeff);
	double *rvector = NULL;

	/* Init the v array */
	unsigned int **v;
	/* v = (unsigned int **) R_alloc(sobolMaxDim, sizeof(unsigned int)); */
	/* for(i=0; i<sobolMaxDim; i++) v[i] = (unsigned int *) R_alloc(sobolMaxCol, sizeof(unsigned int)); */
	v = (unsigned int **) malloc(sobolMaxDim * sizeof(unsigned int));
	for(i=0; i<sobolMaxDim; i++) v[i] = (unsigned int *) malloc(sobolMaxCol * sizeof(unsigned int));
        lastpoint = (unsigned int *) malloc (d * sizeof(unsigned int)); /* (unsigned int *) R_alloc(d, sizeof(unsigned int)); */

        memset (lastpoint, 0, sizeof(unsigned int) * d);
        init_array_v(v, n, d);

	/* Generate a shift */
        if(randomize){
		rvector = (double *) malloc (d * sizeof(double)); /* (double *) R_alloc(d, sizeof(double)); */
		/* GetRNGstate(); */
		for (i = 0; i < d; i++) { *(rvector+i) = 1.0; } /* unif_rand(); } */
		/* PutRNGstate(); */
	}

        /* Compute the recipd */
	numcols = ceil(log(n)/log(2));
        maxn =  ((unsigned int)(1 << numcols));
        recipd = ((double) ((unsigned int)(1 << numcols)));
        recipd = ((double) (1.0 / recipd));

	/* Init result (and randomization) */
        for(i=0; i<d; i++){
		res[i*n] = 0.0;
		if(randomize){
			point = *(lastpoint+i);
			U = *(rvector+i);
			U *= rmaxint;
			randint = (unsigned int) U;
			point = point ^ randint;
			res[i*n] = ((double) point)/rmaxint;
		}
	}

	/* Main loop */
        for(count=0; count<n-1; count++){
		column=0;
		while ((count >> column) & 0x1) /* '>>' = bitwise right shift by column-many bits and fill with 0s from the left; hexadecimal value of 1 */
			column++;
		/* Calculate the next point using Gray Code */
		for (i = 0; i < d; i++)
		{
			*(lastpoint+i) = *(lastpoint+i) ^ v[i][column];
			if(randomize){
				point = *(lastpoint+i);
				/* point = point * (rmaxint*recipd); */
				point = point << (rmaxcoeff - numcols);
				U = *(rvector+i);
				U *= rmaxint;
				randint = (unsigned int) U;
				point = point ^ randint;
				res[i*n+count+1] = ((double) point)/rmaxint;
			}
			else res[i*n+count+1] = ((double) *(lastpoint+i)) * recipd;
		}
	}

	free(v);
	free(lastpoint);
	free(rvector);
}

/* /\** */
/*  * @title R Interface to C for Generating a Sobol Sequence */
/*  * @param n number of points */
/*  * @param d dimension */
/*  * @param randomize logical indicating whether a digital shift should be */
/*           included */
/*  * @return (n, d)-matrix */
/*  * @author Marius Hofert */
/*  *\/ */
/* SEXP sobol_(SEXP n, SEXP d, SEXP randomize) */
/* { */
/* 	/\* Input parameters *\/ */
/* 	int n_ = asInteger(n); /\* numeric(1) *\/ */
/* 	int d_ = asInteger(d); /\* numeric(1) *\/ */
/* 	int randomize_ = asLogical(randomize); /\* 0 (= FALSE), 1 (= TRUE) *\/ */

/* 	/\* Create result object *\/ */
/* 	SEXP res = PROTECT(allocMatrix(REALSXP, n_, d_)); /\* (n,d)-matrix *\/ */
/* 	double *res_ = REAL(res); /\* pointer to the values of res *\/ */

/* 	/\* Main *\/ */
/* 	sobol(n_, d_, randomize_, res_); */

/* 	/\* Return *\/ */
/* 	UNPROTECT(1); /\* clean-up *\/ */
/* 	return res; */
/* } */


int main() {
	int n=10; /* 10000000;*/
	int d=5;
	int randomize=0;
	double *res = (double *) malloc(n*d * sizeof(double));
	sobol(n, d, randomize, res);
	/* sobol(n, d, randomize, res); */
	/* int k; */
	/* for(k=0; k<n*d; k++) { */
	/* 	printf("%f\n", res[k]); */
	/* } */

	int i,j;
	for(j=0; j<d; j++) {
		for(i=0; i<n; i++) {
			printf("%f ", res[i+j*n]);
		}
		printf("\n");
	}
}