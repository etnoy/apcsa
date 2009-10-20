#define LOOPS 1000
#define BUFSIZE 100000
#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main (int argc, char** argv) {
	int rank, size,i,t;

	MPI_Status stat;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double sendbuf[BUFSIZE];
	double recvbuf[BUFSIZE];
	int buf;
	double timestamp;
	double bandwidth;
	int msglen;

	for(int j=0;j<19;j++) {
		buf=pow(2,j);
		timestamp=MPI_Wtime();
		for(int i=0;i<LOOPS; i++) {
			if(rank==0)
			{
				MPI_Send(&sendbuf, buf, MPI_DOUBLE, size-1,1,MPI_COMM_WORLD);
				MPI_Recv(&sendbuf, buf, MPI_DOUBLE, size-1,1,MPI_COMM_WORLD,&stat);
			} else if (rank==size-1) {
				MPI_Recv(&recvbuf, buf, MPI_DOUBLE, 0,1,MPI_COMM_WORLD,&stat);
				MPI_Send(&recvbuf, buf, MPI_DOUBLE, 0,1,MPI_COMM_WORLD);
			}
		}
		if(rank ==0)
		{
			timestamp=MPI_Wtime()-timestamp;
			timestamp=MPI_Wtick()*timestamp;
			bandwidth=buf*LOOPS/timestamp/2;
			printf("Size: %d: %f Bps\n", buf, bandwidth/8);
		}
	}


	MPI_Finalize();

	return 0;

}
