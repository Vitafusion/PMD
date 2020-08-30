#include <R.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <tgmath.h>
#include <complex.h>
#include <fftw3.h>


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








