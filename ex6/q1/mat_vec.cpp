#include <iostream>
#include "walltime.h"
#include <mpi.h>
#define N 10003					  // Each processor will not have equal number of rows. Try to make the distribution as equal as possible

using namespace std;

int main (int argc, char** argv)
{

	int rank,size;

	// Initialize MPI 
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int rows=0, start_row=0;			// rows is the number of rows assigned to each processor and start_row is the global row index of the first row on each processor
	double *u = new double[N];		// u is the L.H.S vector that is replicated on all processors
	double *A, *f;                    // A is NxN only on rank 0 and the size of A and f on other processors depends on the distribution
	double* f_gather;		        	// This is the resultant vector allocated only on processor 0; 	
	double time=0.0,start_time=0.0;
	int *displ, *recv_count;			// You will need these arrays for MPI_Gatherv(...)
	int stepSize = (N-N%(size-1))/(size-1);
	int rows_buf[1];
	int start_row_buf[1];

	if (rank == 0)
	{
		// 1.a) Allocate memory for matrix and generate the hilbert matrix
		A = new double[N*N];	 	
		for (int i = 0; i < N; i++)
			for (int j=0; j <N;j++)
				A[i*N+j] = 1.0/(i+j+1.0);

		// 1.b) Generate L.H.S vector with all elements as 1.0
		for(int i=0; i < N; i++)
			u[i] = 1.0;

		// 2.a) Calculate the rows to be distributed to each processors and distribute the matrix and other relevant information
		displ = new int[size];
		recv_count = new int [size];

		displ[0] = 0;
		recv_count[0] = 0;

		for(int jj = 1; jj<size; jj++){
			start_row = (jj-1)*stepSize;

			if(rank != size - 1)
				rows = stepSize;
			else
				rows = N-start_row;

			rows_buf[0] = rows;
			start_row_buf[0] = start_row;

			MPI_Send(&rows_buf, 1, MPI_INT, jj, jj+1, MPI_COMM_WORLD);	
			MPI_Send(&start_row_buf, 1, MPI_INT, jj, jj+2, MPI_COMM_WORLD);	
			MPI_Send(&A[start_row*N], rows*N, MPI_DOUBLE, jj, jj, MPI_COMM_WORLD);
			
			displ[jj] = start_row;
			recv_count[jj] = rows;
		}

		f = new double[rows];	
		f_gather = new double [N];
	}
	else
	{   
		//2.b) Receive row-blocks and other relevant information from master processor
		MPI_Recv(&rows_buf, 1, MPI_INT, 0, rank+1, MPI_COMM_WORLD, &status);	
		MPI_Recv(&start_row_buf, 1, MPI_INT, 0, rank+2, MPI_COMM_WORLD, &status);	
		
		rows = rows_buf[0];
		start_row = start_row_buf[0];

		A = new double[rows*N];

		MPI_Recv(A, rows*N, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);

		f = new double[rows];
	}   

	//3) Replicate L.H.S matrix on each processor

	MPI_Bcast(u, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

	cout << "My rank = " << rank << ", rows = " << rows << ", start row = " << start_row << endl;  	

	// 4) Calculate part of the resultant vector

	double temp;
	if(rank!=0)
		for(int jj = 0; jj<rows; jj++){
			temp = 0.0;
			for(int ii = 0; ii<N; ii++)
				temp += A[ii+jj*N]*u[ii];
			f[jj] = temp;
		}

	//5) Gather the resultant vector on the master processor 

	MPI_Gatherv(f,rows,MPI_DOUBLE,f_gather,&recv_count[rank],&displ[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(rank == 0)
	{
		cout << "Execution time for problem size " << N*N << " using " << size << " processors is "  << time << endl;
	}

	// Delete dynamically allocated memory
	delete [] A;
	delete [] u;
	delete [] f;

	if(rank == 0)
	{
		delete [] f_gather;
		delete [] recv_count;
		delete [] displ;
	}

	MPI_Finalize(); 
	return 0;

}
