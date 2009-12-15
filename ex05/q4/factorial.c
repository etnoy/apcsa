#define BUFSIZE 100
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main (int argc, char** argv) {
	int rank, size,i,t,sum;
	MPI_Status stat;
	MPI_Request sendreq, recvreq;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int sendbuf[1];
	int recvbuf[1];

	sendbuf[0]=rank+1;	
	MPI_Scan(&sendbuf, &recvbuf, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
	sum=recvbuf[0];
	printf("Result: %d\n", sum);

}
