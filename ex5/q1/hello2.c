#define bufsize 100
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main (int argc, char** argv) {
	int rank, size,i,t;
	char greeting[bufsize];

	MPI_Status stat;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	if(rank==0)
	{
		sprintf(greeting,"Hello world from process %d of %d\n", rank,size);
		fputs(greeting,stdout);
		for(int partner=1;partner<size; partner++) {
			MPI_Recv(greeting, sizeof(greeting), MPI_BYTE, partner,1,MPI_COMM_WORLD,&stat);
			fputs(greeting,stdout);
		}
	}else {
		printf("Process %d of %d\n", rank,size);
		MPI_Send(greeting, strlen(greeting)+1,MPI_BYTE,0,1,MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;

}
