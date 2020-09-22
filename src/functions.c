#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <tgmath.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>
#define MAXCOL 10
#define BINS 1


//*****************************************************************************//
void l_vec_compute(int k, int *l_vec, int *cn_vec, int m)
{
  int i, aa, bb;

  for(i=0; i<m-1; i++)
  {
    aa=k%cn_vec[i];
    bb=(k-aa)/cn_vec[i];
    l_vec[i]=bb;
    k=aa;
  }

  return;
}


//*****************************************************************************//
void pmn_mdfft(double *res, int *nnt, int *nn, int *mm, double *pp, int *nn_vec, int *l_vec, int *cn_vec)
{
  fftw_complex *in, *out;
  int i, j, k;
  int n, m, nt;
  fftw_plan p;
  double tmp, con, pij, pim, ww;
  double complex ctmp, ctmp1, ctmp2, qval, a1, a2;

  nt=nnt[0]; //(n+1)^(m-1)
  n=nn[0]; // draws
  m=mm[0]; // categories

  //Rprintf("nt %u, n %u, m %u \n", nt, n, m);

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nt);

  ww=2*PI/(n+1);

  //Rprintf("ww %lf \n", ww);

  for(k=0; k<nt; k++)
  {
    qval=0.0 + 0.0 * I;

    l_vec_compute(k, l_vec, cn_vec, m);

    //for(ii=0; ii<m-1; ii++)
    //{
    //  Rprintf("l_vec %u \n", l_vec[ii]);
    //}

    for(i=0; i<n; i++)
    {
      ctmp=0.0 + 0.0 * I;

      for(j=0; j<m-1; j++)
      {
        pij=pp[i+n*j];

        //printf("pij: %lf, l_vec[j], %u, ww, %lf, \n", pij, l_vec[j], ww);

        ctmp1=0.0+l_vec[j]*ww*I;

        //printf("ctmp1: %lf +%lf*i\n", creal(ctmp1), cimag(ctmp1));

        a1=pij+0.0*I;
        a2=exp(ctmp1);
        //printf("a1: %lf +%lf*i\n", creal(a1), cimag(a1));
        //printf("a2: %lf +%lf*i\n", creal(a2), cimag(a2));

        ctmp2=a1*a2;

        //printf("ctmp2: %lf +%lf*i\n\n", creal(ctmp2), cimag(ctmp2));

        ctmp+=ctmp2;
      }

      pim=pp[i+n*(m-1)];
      ctmp+=pim;

      ctmp=log(ctmp);
      qval+=ctmp;
    }

    qval=exp(qval);

    in[k]=qval;

    //printf("qval: %lf +%lf*i\n", creal(qval), cimag(qval));
    //printf("in[k]: %lf +%lf*i\n", creal(in[k]), cimag(in[k]));



  }

  p=fftw_plan_dft(m-1, nn_vec, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);

  con=pow(n+1, m-1);

  for(k=0; k<nt; k++)
  {
    tmp=creal(out[k]);
    res[k]=tmp/con;
    //printf("out[k]: %lf +%lf*i\n", creal(out[k]), cimag(out[k]));
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return;
}
//*****************************************************************************//


void pmd_simulation(double *res0, int *nnt, int *nn, int *mm, double *pp, int *nn_vec, int *l_vec, int *cn_vec, double d)
{
  int nt,i,j,k;
  const unsigned int sum_n = BINS;
  const unsigned int m=mm[0];
  nt=nnt[0]; //(n+1)^(m-1)
  int (*sim)[m] = malloc(d * sizeof * sim);
  double (*p)[m] = malloc(nn[0] * sizeof(double) * m);
  for(i=0;i<nn[0];i++)
	for(j=0;j<m;j++)
		p[i][j]=pp[i+nn[0]*j];

  //for(i=0;i<nn[0];i++)
	//for(j=0;j<m;j++)
   	//	printf("%lf",p[i][j]);


    /*gsl variable */
  const gsl_rng_type * T;
  gsl_rng * r_global;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r_global = gsl_rng_alloc (T);
  int sum[m];
  unsigned int n[m];
  for(i=0;i<m;i++){
	sum[i]=0;
  }
/*  for(i=0;i<m;i++){
	printf("%d\n",sum[i]);
  }
*/

/* creat simulation matrix */
  for(k=0;k<d;k++){
  	for(i=0;i<m;i++){
		sum[i]=0;
  	}
	for(i=0;i<nn[0];i++){
		gsl_ran_multinomial ( r_global, m, sum_n, p[i], n);
		for(j=0;j<m;j++){
			//printf("%d\n",n[j]);
        		sum[j]+=n[j];
		}

        }
	//for(j=0;j<m;j++)
	//	printf("%d\n",sum[j]);
	for(i=0;i<m;i++){
		sim[k][i]=sum[i];
	}
  }

  //for(i=0;i<d;i++){
//	for(j=0;j<m;j++){
//		printf("%d\n",sim[i][j]);
//	}
 // }

//compare to giving l_vec,x_vex
  int x_vec[m];
  int u,v,temp;
  double count;
  int flag;
  for(i=0;i<nt;i++){
  	l_vec_compute(i, l_vec, cn_vec, m);
	temp=0;
	count=0;
	//for(int ii=0; ii<m-1; ii++)
    	//{
      	//	printf("l_vec %u \n", l_vec[ii]);
    	//}
	for(j=0;j<m-1;j++){
		x_vec[j]=l_vec[j];
        }
        for(k=0;k<(m-1);k++){
		temp+=l_vec[k];
	}

	if(temp < nn[0]){
		x_vec[m-1]=nn[0]-temp;
	}
	else{
		x_vec[m-1]=0;
	}
        //printf("temp %d\n",temp);

	//for(j=0;j<m;j++){
	//	printf("x_vec: %u \n",x_vec[j]);
        //}
	for(u=0;u<d;u++){
		flag=1;
		for(v=0;v<m;v++){
			if(sim[u][v]!=x_vec[v])
				flag=0;
		}
		if(flag==1)
			count++;
	}
	res0[i]= count / d;
  }
 /* for(i=0;i<nt;i++){
	printf("%lf\n",res0[i]);
  }
*/
return;
}
//*****************************************************************************//

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



