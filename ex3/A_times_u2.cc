#include <omp.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "crs.h"
#include "walltime.h"

extern "C" {
    double ddot_(const int* n, const double* dx, const int* incx, const double* dy, const int* incy);
}

void reset_vector(double vector[], int dim);
void equality_test(int N, double v[], double v_test[]);

void A_dense_mult_ordinary(int N, double A[], double u[], double v[]);
void A_dense_mult_chunk   (int N, double A[], double u[], double v[]);
void A_dense_mult_BLAS    (int N, double A[], double u[], double v[]);

void A_triangular_mult_ordinary (int N, double A[], double u[], double v[]);
void A_triangular_mult_optimized(int N, double A[], double u[], double v[]);
void A_triangular_mult_self_opt (int N, double A[], double u[], double v[]);

using namespace std; 

int main(int argc, char** argv){
    int    p, row, column, N = 4096, N_sparse = 100000;
    double time, time_start = 0.0;

    double * u         = new double[N];       assert( u         != 0 );
    double * u_sparse  = new double[N_sparse];assert( u_sparse  != 0 );
    double * v1        = new double[N];       assert( v1        != 0 );
    double * v2        = new double[N_sparse];assert( v2        != 0 );
    double * v_dense   = new double[N];       assert( v_dense   != 0 );
    double * v_triag   = new double[N];       assert( v_triag   != 0 );
    double * v_sparse  = new double[N_sparse];assert( v_sparse  != 0 );
    double * v_sparseT = new double[N_sparse];assert( v_sparseT != 0 );
    double * A         = new double[N*N];     assert( A         != 0 );
    CRSmat * A_sparse  = new CRSmat();

    reset_vector(v1       , N);        // v1 : for parallel calculations (dense)
    reset_vector(v2       , N_sparse); // v2 : for parallel calculations (sparse)
    reset_vector(v_dense  , N);
    reset_vector(v_triag  , N); 
    reset_vector(v_sparse , N);
    reset_vector(v_sparseT, N);

    for (int i=0; i<N; i++)
        u[i] = ((double)i);

    for (int ii=0; ii<N_sparse; ii++)
        u_sparse[ii] = ((double)ii);

    for (int iii=0; iii<N*N; iii++)
        A[iii] = ((double)iii); 

    // generating a sparse matrix
    A_sparse->beginAssembly(N_sparse, N_sparse);
    double fillfactor = 0.05;
    int    gap        = (int)(N_sparse * fillfactor);
    for (row=0; row<N_sparse; row++){
        for (column=0; column<N_sparse; column += gap){
            A_sparse->setAij(row, column, fillfactor*column);
        }
    }


    #pragma omp parallel shared(p)
    {
        p = omp_get_num_threads();
    }
    cout << endl << "Number of threads: " << p << endl << endl;


    // a)

    cout << "a) dense A           (Serial computation of A_dense*u = v_dense: ";

    // calculate test vector: v_dense (serial multiplication)
    time = walltime(&time_start);
    for (row=0; row<N; row++){
        double sum = 0.0;
        for (column=0; column<N; column++)
            sum += A[row*N+column]*u[column];
        v_dense[row] = sum;
    }
    time = walltime(&time);
    cout << time << " sec )" << endl << endl;


    cout << "   > ordinary  : ";
    time = walltime(&time_start);
    A_dense_mult_ordinary(N, A, u, v1);        
    time = walltime(&time);
    equality_test(N, v1, v_dense);
    cout << time << " sec " << endl;

    cout << "   > stat chunk: ";
    reset_vector(v1, N);
    time = walltime(&time_start);
    A_dense_mult_chunk(N, A, u, v1);           // to implement
    time = walltime(&time); 
    equality_test(N, v1, v_dense);
    cout << time << " sec " << endl;

    cout << "   > BLAS      : ";
    reset_vector(v1, N);
    time = walltime(&time_start);
    A_dense_mult_BLAS(N, A, u, v1);            // to implement
    time = walltime(&time);
    equality_test(N, v1, v_dense);
    cout << time << " sec " << endl;



    // b)

    cout << endl << "b) triangular A      (Serial computation of A_dense*u = v_triag: ";

    // calculate test vector: v_triag (serial multiplication)
    time = walltime(&time_start);
    for (row=0; row<N; row++){
        double sum = 0.0;
        for (column=row; column<N; column++)
            sum += A[row*N+column]*u[column];
        v_triag[row] = sum;
    }
    time = walltime(&time);
    cout << time << " sec )" << endl << endl;


    cout << "   > ordinary  : ";
    reset_vector(v1, N);
    time = walltime(&time_start);
    A_triangular_mult_ordinary(N, A, u, v1);   // to implement
    time = walltime(&time);
    equality_test(N, v1, v_triag);
    cout << time << " sec " << endl;

    cout << "   > dyn. chunk: ";
    reset_vector(v1, N);
    time = walltime(&time_start);
    A_triangular_mult_optimized(N, A, u, v1);   // to implement
    time = walltime(&time);
    equality_test(N, v1, v_triag);
    cout << time << " sec " << endl;

    cout << "   > self opt. : ";
    reset_vector(v1, N);
    time = walltime(&time_start);
    A_triangular_mult_self_opt(N, A, u, v1);    // to implement
    time = walltime(&time);
    equality_test(N, v1, v_triag);
    cout << time << " sec " << endl;




    // c)
    cout << endl << "c) sparse A          ";
    A_sparse->completeAssembly();
    cout << "   > serial    :                     ";
    time = walltime(&time_start);
    A_sparse->mult_serial(u_sparse, v_sparse);  // to implement in "crs.h"
    time = walltime(&time);
    cout << time << " sec " << endl;


    cout << "   > parallel  : ";
    time = walltime(&time_start);
    A_sparse->mult_OMP(u_sparse, v2);           // to implement in "crs.h"
    time = walltime(&time);
    equality_test(N_sparse, v2, v_sparse);
    cout << time << " sec " << endl;





    // d)

    cout << endl << "d) sparse transposed A" << endl << endl;
    cout << "   > parallel  :                     ";
    reset_vector(v2, N_sparse);
    time = walltime(&time_start);
    A_sparse->mult_trans_OMP(u_sparse, v2);    // to implement in "crs.h"
    time = walltime(&time);
    cout << time << " sec " << endl;


    return 0;

}


void reset_vector(double vector[], int dim){
    int i;
#pragma omp for  
    for (i=0; i<dim; i++)
        vector[i] = 0.0;
}

void equality_test(int N, double v[], double v_test[]){
    cout.width(20);  
    for (int i=0; i<N; i++){

        if (v[i] != v_test[i]){
            cout << left << "v != v_test !!!!!!";
            return;
        }
    }
    cout << left << "v = v_test";
}


// ------------- IMPLEMENT YOUR ROUTINES HERE ------------------------


// a)
void A_dense_mult_ordinary(int N, double A[], double u[], double v[]){
	
	#pragma omp parallel for
	
	for(int i = 0; i<N; i++)
	{
		double temp_sum = 0.0;
		for(int j = 0; j<N; j++)
			temp_sum += A[ i*N + j ]*u[j];
		v[i] = temp_sum;
	}	

	return;
}


void A_dense_mult_chunk(int N, double A[], double u[], double v[]){
	
	// Here we edit to get chunks of size half the equal, equal and double the equal:
	int chunk = N/omp_get_num_threads(); 

	#pragma omp parallel for schedule(static, chunk)
	
	for(int i = 0; i<N; i++)
	{
		double temp_sum = 0.0;
		for(int j = 0; j<N; j++)
			temp_sum += A[ i*N + j ]*u[j];
		v[i] = temp_sum;
	}	

	return;
}



void A_dense_mult_BLAS(int N, double A[], double u[], double v[]){

 	const int I = 1;	
	#pragma omp parallel for 
		for(int i = 0; i<N; i++)
		{
			v[i] = ddot_(&N,&A[i*N],&I,u,&I);
		}	
	return;
}



// b)

void A_triangular_mult_ordinary(int N, double A[], double u[], double v[]){
	
	double temp_sum = 0.0;
	
	#pragma omp parallel for private(temp_sum)
	for(int i = 0; i<N; i++)
	{
		for(int j = i; j<N; j++)
			temp_sum += A[ i*N + j ]*u[j];
		v[i] = temp_sum;
	}	
	return;


}


void A_triangular_mult_optimized(int N, double A[], double u[], double v[]){

	double temp_sum = 0.0;
	// Here we edit to get chunks of size half the equal, equal and double the equal:
	int chunk = ( int ) (ceil(N*(N+1)/2/omp_get_num_threads())); 

	#pragma omp parallel for private(temp_sum) schedule(dynamic, chunk)
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<N; j++)
			temp_sum += A[ i*N + j ]*u[j];
		v[i] = temp_sum;
	}	
	return;

}


void A_triangular_mult_self_opt (int N, double A[], double u[], double v[]){



}

