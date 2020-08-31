
/*
	gcc -Wall -Werror -fPIC -I/home/admin1/gsl/include -I"/usr/share/R/include" -c simulation.c	
gcc -shared -L/home/admin1/gsl/lib -o simulation.so simulation.o -lgsl -lgslcblas -lm
*/
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include<math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define MAXCOL 10
#define BINS 1

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
