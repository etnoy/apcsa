#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "walltime.h"

#define N_R      1000
#define N_I      1000
#define V_R_MIN -1.0
#define V_R_MAX  1.0
#define V_I_MIN -1.0
#define V_I_MAX  1.0
#define CMAX 10000

static int mandelbrot(float, float);

int main()
{
  int     field[N_R*N_I];
  int     dither[N_R*N_I];
  float   v_r, v_i;
  float   v_r_step, v_i_step;
  int     i,j,myId;
  int     num_threads;
  double time, time_start=0.0;
  int counter;
  int N;
  N=omp_get_thread_num();
  size=omp_get_num_threads();

  v_r_step = (V_R_MAX-V_R_MIN)/N_R;
  v_i_step = (V_I_MAX-V_I_MIN)/N_I;

  time = walltime(&time_start);
  for( j=1; j<=N_I; j++) {                               // iterate over the columns
      v_i = V_I_MIN+(j-1)*v_i_step;
      for( i=1; i<=N_R; i++ ) {                          //iterate over the rows
          v_r = V_R_MIN+(i-1)*v_r_step;
          field[i-1 +(j-1)*N_R] = mandelbrot(v_r, v_i);  // column major
      }
    }

  for( j=0; j<N_I; j++) 
     for( i=1; i<N_R-1; i++ )
 	dither[i+j*N_R] = 0.5*field[i+j*N_R] + 0.25*(field[(i-1)+j*N_R] + field[(i+1)+j*N_R]); 


  time = walltime(&time);
  printf("Serial execution time = %f sec\n",time);



  // Use can use the following piece of code to verity the fractal pattern for small image sizes 

  for( i=0; i<N_R; i++)
  {  for ( j=0; j<N_R; j++)
	if(field[i+j*N_R] == -1)
	   printf("*"); 
        else 
           printf(" ");
        printf("\n");
  }

  return 0;
}

static int mandelbrot(float vr, float vi)
{
  struct { float r,i; } cc, zz, tmp;
  int i;

  cc.r = vr;
  cc.i = vi;
  zz.r = 0.0;
  zz.i = 0.0;
  for( i=1; i<=CMAX; i++ ) {
    tmp.r = (zz.r*zz.r-zz.i*zz.i)+cc.r;
    tmp.i = (zz.r*zz.i+zz.i*zz.r)+cc.i;
    zz.r = tmp.r;
    zz.i = tmp.i;
    if((zz.r*zz.r)+(zz.i*zz.i) > 5.0 )
      return i;
  }
  return -1;
}
