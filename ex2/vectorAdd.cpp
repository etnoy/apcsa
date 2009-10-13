#include <omp.h>
#include <iostream>
#include "walltime.h"
#define NMAX 1000000

using namespace std; 

int main()
{	int myId, numTdreads;
	double time, time_start=0.0;
	double dotProduct;
	double *a,*b,*c,d;
	int chunk_size, n_start,n_end;

	// Allocate memory for the vectors as 1-D arrays
	a = new double[NMAX];
	b = new double[NMAX];
	c = new double[NMAX];

	
	// Initialize the vectors with some values
	for(int i=0; i<NMAX; i++)
	{
		a[i] = i;
		b[i] = i/10.0;
	}

	// serial execution
	time = walltime(&time_start);
	for(int i=0; i< NMAX; i ++)
	{
		c[i] = 2*a[i] + b[i];
	}
	time = walltime(&time);
	cout << "Serial execution time =" << time << "sec" << endl;

	time = walltime(&time_start);

	// parallel execution
	
	// TODO: Modify the following code to obtain sppedup (minor change required)
	omp_set_num_threads(4);
	chunk_size=NMAX/omp_get_num_threads();
	cout << omp_get_num_threads() << endl;

	#pragma omp parallel private(time)
	{
	for(int i=omp_get_thread_num()*chunk_size; i< (1+omp_get_thread_num())*chunk_size && i< NMAX; i++)
    {
	  //  printf("Thread %d: i=%d\n", omp_get_thread_num(),i);
    	c[i] = 2*a[i] + b[i];
    }
	}
    time = walltime(&time);
    cout << "Parallel execution time =" << time << "sec" << endl;

	// TODO: Implement scalar procduct (dot product)

	#pragma omp parallel private(time)
	{
	for(int i=omp_get_thread_num()*chunk_size; i< (1+omp_get_thread_num())*chunk_size && i< NMAX; i++)
    {
	  //  printf("Thread %d: i=%d\n", omp_get_thread_num(),i);
    	d = d +a[i]*b[i];
    }
	}

    time = walltime(&time);
    cout << "Parallel execution time =" << time << "sec" << endl;
	// De-allocate memory
	delete [] a;
	delete [] b;
	delete [] c;	
	
	return 0;
}
