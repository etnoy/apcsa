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
	int receiver,sender;

	if (rank==size-1 )
	{
		receiver=0;
	} else {
		receiver=rank+1;
	}

	if (rank==0)
	{
		sender=size-1;
	} else {
		sender=rank-1;
	}

	for(int i=0;i<size;i++) {
		sendbuf[0]=rank;	
		MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		sum=recvbuf[0];
	}
	printf("Sum: %d\n", sum);

}
