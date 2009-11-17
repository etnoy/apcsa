#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<string.h>

unsigned int MATRIX_SIZE = 100; //for the beginning

typedef struct {
	int rank; // Rank of the communicator in COMM_WORLD
	int size; // How many communicators in COMM_WORLD
	int q; // nb of submatrices per row/column.
	int row_index; // the sub_matrix row-index
	int col_index; // the sub_matrix col-index
	int row_start; // the index of the first row
	int col_start; // the index of the first col
	MPI_Comm col_comm; // the column communicator (send upwards)
	MPI_Comm row_comm; // the row communicator (broadcast).
	int col_comm_rank; // rank in the column communicator
	int row_comm_rank; // rank in the row communicator
	double* local_A; // this is calculated in the beginning and remains const
	double* local_B; // this is calculated in the beginning and remains const 
	// local_B is actually not really necessary.
	double* A; 
	double* B;
	double* C;
} COMM_INFO;

// This function fills in all the elements of comm_info (not the matrices)
void generate_communicators(COMM_INFO* comm_info) {
	//

	MPI_Comm_rank(MPI_COMM_WORLD, &comm_info->rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_info->size);

	int *dims;
	int *pers;

	dims = malloc(comm_info->size*sizeof(int));
	pers = malloc(comm_info->size*sizeof(int));
	int newsize=MATRIX_SIZE/comm_info->size;

	for(int i=0;i<comm_info->size;i++)
	{
		dims[i]=newsize;
		pers[i]=1;
	}

	MPI_Cart_create(comm_info->col_comm, MATRIX_SIZE, dims,pers ,0, &(comm_info->col_comm));
	MPI_Cart_create(comm_info->row_comm, MATRIX_SIZE, dims,pers ,0, &(comm_info->row_comm));

	return;

}

// This function allocates and fills in A and B according to the exercise sheet
double* generate_sub_matrix(COMM_INFO* comm_info) {

	int submatr_size;
	int rows;
	double *matrx;


	submatr_size=MATRIX_SIZE/comm_info->size;

	rows=sqrt(comm_info->size);
	matrx= malloc(submatr_size*submatr_size*sizeof(double));

	for(int i=comm_info->row_start;i<comm_info->row_start+submatr_size;i++)
	{
		for(int j=comm_info->col_start;j<comm_info->col_start+submatr_size;j++)
		{
			matrx[(i-comm_info->row_start)*rows+(j-comm_info->col_start)]=(cos(i*rows+j));

		}
	}

	return matrx;

}

// The process that contains the global element (i,j) sends it to the master 
// process. The master will print it to the console.
void collect_and_print_element(COMM_INFO* comm_info,
		unsigned int i, unsigned int j) {

	if(comm_info->rank==0)
	{
		MPI_Irecv( ///todo: jobba här och fixa en recv non-blocking :) //J
	} else {
	}

}

// This function adds up the sub-matrix C.
int multiply_local(COMM_INFO* comm_info) {

	int submatr_size;
	int rows;
	submatr_size=MATRIX_SIZE/comm_info->size;

	rows=sqrt(comm_info->size);
	for(int i=comm_info->row_start;i<comm_info->row_start+submatr_size;i++)
	{
		for(int j=comm_info->col_start;j<comm_info->col_start+submatr_size;j++)
		{
			for(int k=comm_info->col_start;k<comm_info->col_start+submatr_size;k++)
			{
				comm_info->C[i,j]=comm_info->C[i,j]+comm_info->local_A[(i-comm_info->row_start),k]*comm_info->local_B[k,(j-comm_info->col_start)]; 
			}
		}
	}
	return 4; //really random number, guranteed by dice roll

}


int main(int argc, char** argv) {

	COMM_INFO comm_info;
	MPI_Init(&argc, &argv);

	generate_communicators(&comm_info);

	int q = comm_info.q;
	int sub_matrix_size = MATRIX_SIZE / q;

	comm_info.C = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));
	comm_info.A = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));
	comm_info.B = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));

	//local_A and local_B will be allocated in generate_sub_matrix:
	comm_info.local_A = generate_sub_matrix(&comm_info);
	comm_info.local_B = generate_sub_matrix(&comm_info);

	memcpy(comm_info.B, 
			comm_info.local_B, 
			sub_matrix_size * sub_matrix_size * sizeof(double));
	unsigned int i;

	//initialize C to 0: // see also memset
	for(i = 0; i < sub_matrix_size*sub_matrix_size; i++) {
		comm_info.C[i] = 0.0;
	}

	///////////////////////////////////////////////////////////////////
	// TODO: subquestion e)
	///////////////////////////////////////////////////////////////////


	int rank_to_send_to; // ... 
	int rank_to_get_from;

	// ...

	collect_and_print_element(&comm_info,98,99);

	free(comm_info.A);
	free(comm_info.B);
	free(comm_info.C);
	free(comm_info.local_A);
	free(comm_info.local_B);

	MPI_Finalize();

}
