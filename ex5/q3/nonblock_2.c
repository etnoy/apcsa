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
	int ranks[size];

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
	sum=0;

	ranks[0]=rank;	
	for(int i=0;i<size;i++) {
		sendbuf[0]=ranks[i];	
		MPI_Sendrecv(&sendbuf, 1, MPI_INT, receiver,1, &recvbuf, 1, MPI_INT, sender, 1,MPI_COMM_WORLD, &stat);
		ranks[(i+1)%size]=recvbuf[0];
		sum+=recvbuf[0];
	}
	printf("Sum: %d\n", sum);

}
