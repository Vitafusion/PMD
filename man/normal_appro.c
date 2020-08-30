
/*
gcc -Wall -Werror -fPIC -I/home/admin1/gsl/include -I"/usr/share/R/include" -c normal_appro.c
gcc -shared -L/home/admin1/gsl/lib -o normal_appro.so normal_appro.o -lgsl -lgslcblas -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void pmd_normal(int *nn, int *mm, double *pp, int *x_vec)
{
  int n,i,j,k;
  n=nn[0];
  size_t m=mm[0];

//input pp
  double (*p)[m] = malloc(nn[0] * sizeof(double) * m);
  for(i=0;i<nn[0];i++)
	for(j=0;j<m;j++)
		p[i][j]=pp[i+nn[0]*j];

//define result variable
  double obs_res = 0;

//gsl variables
  gsl_vector * mu = gsl_vector_calloc(m);
  gsl_vector * temp0 = gsl_vector_calloc(m);
  gsl_matrix * Sigma = gsl_matrix_calloc(m,m);
  gsl_matrix * temp = gsl_matrix_calloc(m,m);
  gsl_matrix * L = gsl_matrix_calloc(m, m);
  gsl_vector * x = gsl_vector_calloc(m);
  gsl_vector * work = gsl_vector_calloc(m);

//initialize x as values in x_vec

  for(i=0;i<m;i++){
	gsl_vector_set(x,i,x_vec[i]);
  }
// show x
//  for(i=0;i<m;i++){
//	printf("x_vec[%d] = %g\n", i, gsl_vector_get (x, i));
//  }
//
  for(i=0;i<m;i++){
	for(j=0;j<m;j++){
		gsl_matrix_set(Sigma,i,j,0);	
	}
  }
  for(i=0;i<m;i++){
	gsl_vector_set(mu,i,0);
  }

//calculate covariance matrix;
  for(k=0;k<n;k++){
  	for(i=0;i<m;i++){
		for(j=0;j<m;j++){
			//printf("%lf\n",p[0][j]);
			//printf("%lf\n",p[0][i]);
			if(i==j)
				gsl_matrix_set(temp,i,j,(p[k][i]-p[k][i]*p[k][j]));
				
				
			else
				gsl_matrix_set(temp,i,j,(0-p[k][i]*p[k][j]));
	
		}
  	}
	gsl_matrix_add(Sigma,temp);
  }
//calculate mu vector
  for(k=0;k<n;k++){
	for(i=0;i<m;i++){
		gsl_vector_set(temp0,i,p[k][i]);	
	}
	gsl_vector_add(mu,temp0);
  }
//print mu vector
//  for(i=0;i<m;i++){
//	printf("mu_vec[%d] = %g\n", i, gsl_vector_get (mu, i));
//  }
//print sigma matrix
//  for(i=0;i<m;i++)
//	for(j=0;j<m;j++)
//		printf ("m(%d,%d) = %g\n", i, j, gsl_matrix_get (Sigma, i, j));


  gsl_matrix_memcpy(L, Sigma);
  gsl_linalg_cholesky_decomp1(L);

// calculate probability of given x
  gsl_ran_multivariate_gaussian_pdf(x, mu, L, &obs_res, work);
  printf("obs_res: %lf\n",obs_res);

//result:  printf("obs_res: %lf\n",obs_res);
  free(*p);
  gsl_vector_free(mu);
  gsl_matrix_free(Sigma);
  gsl_matrix_free(L);
  gsl_vector_free(x);
  gsl_vector_free(work);
  
  return;
}
