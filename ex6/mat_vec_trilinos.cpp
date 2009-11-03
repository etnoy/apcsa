#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include <iostream>
#define N 10

int main(int argc, char *argv[])
{	
    MPI_Init(&argc,&argv);

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    Epetra_Time Time(Comm);
    int rank = Comm.MyPID();

    int NumGlobalElements = N;
	
    // The following map will distribute the matrix row-wise to different processors
    Epetra_Map map(N,						// Total number of rows
				   0,                       // Array indexing begins from 
				   Comm);					// Communicator

    int numMyElements = map.NumMyElements();                    // Number of rows assigned to the processoe
    int* myGlobalElements = map.MyGlobalElements();             // Row numbers on the processor
   
    Epetra_Vector u(map);					// L.H.S. vector 
    u.PutScalar(1.0);

    Epetra_Vector f(map);					// R.H.S. vector

    Epetra_CrsMatrix A(Copy,				
                       map,					// map specifying the matrix block distribution
                       N);			        // Number of elements in each row
    
    // Note that even though the Hilbert Matrix is dense, it is stored as a CRS matrix, hence we need to specify the column indices in each row
    double *values = new double[N];			//row buffer to be copied into the matrix
    int *indices = new int[N];              //column indices

    for(int i=0; i<N; i++)
	indices[i] = i;

	// Generate the Hilbert matrix
    for(int i=0; i< numMyElements; i++)
    {
       for (int j=0; j<N; j++)
		 values[j] = 1.0/(myGlobalElements[i]+j+1);
       A.InsertGlobalValues(myGlobalElements[i],    // Global row number
 							N, 						// Number or entries in each row
							values,					// values 
							indices);  				// column indices in each row
    }
    
    A.FillComplete();								// This method has to be called to before performing any further operations on the matrix

    A.Multiply(false,         u,     f );
    //         no transpose    L.H.S   R.H.S        //  f = A*u ;  A, u, and f should have the same map 
    cout << f << endl;

    delete [] values;
    delete [] indices;

    MPI_Finalize() ;
    return 0;
}

