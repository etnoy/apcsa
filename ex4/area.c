/*
   This program illustrates the use of Monte Carlo method for calculating the total area 
   of intersecting circles enclosed in a square of unit length.
   Calculating the area analytically is possible but very it is cumbersome
   In this program N points are generated randomly within the bounding square; Then : 
   (Area of circles) / (Area of bounding square) 
   = (No. of points falling inside the circles) / (Total number of points)
   */
//#include <fail.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Define parameters for the circles.
#define X1  0.3
#define Y1  0.5
#define R1  0.3

#define X2  0.55
#define Y2  0.5
#define R2  0.2

#define X3  0.5
#define Y3  0.2
#define R3  0.15

double ggl(double *ds) {
	/* generate u(0,1) distributed random numbers.
	   Seed ds must be saved between calls. ggl is
	   essentially the same as the IMSL routine RNUM.

	   W. Petersen and M. Troyer, 24 Oct. 2002, ETHZ */

	double t,d2=0.2147483647e10;
	t   = *ds;
	t   = fmod(0.16807e5*t,d2);
	*ds = t;
	return((t-1.0e0)/(d2-1.0e0));
}

inline int checkPoint(double x, double y)
{
	//
	// return 1 if point lies within the circumference of any one of the circles
	// else return 0
	if(	sqrt((x-X1)*(x-X1) + (y-Y1)*(y-Y1))<R1 || 
		sqrt((x-X2)*(x-X2) + (y-Y2)*(y-Y2))<R2 || 
		sqrt((x-X3)*(x-X3) + (y-Y3)*(y-Y3))<R3 )
			return 1;
	return 0;
}

int main(int argc, char** argv)
{ 
	// The number of random points is supplied as an argument to main
	if (argc < 2)
	{  printf("Argument missing : Number of random points\n");
		return 0;
	}
	long int N = atoi(argv[1]);
	long int i,count;
	double x,y;
	double ts,te;
	double seed_x,seed_y;
	seed_x = 331.0; seed_y = 151.0;
	double area;  

	count=0;
	ts=omp_get_wtime();

#pragma omp parallel for schedule(runtime)
	for (i=1;i<=N;i++)
	{ 
		x = ggl(&seed_x);
		y = ggl(&seed_y);  
		if(checkPoint(x,y) == 1) 
			count++;
	}

	area = (double)count/N;
	te = omp_get_wtime();

	printf( "Computed area = %4.6f  in  %5.2f seconds \n", area, te-ts );

	return 0;
}
