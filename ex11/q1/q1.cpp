#include <iostream>;

using namespace std;


int bitonic_split_sequential(int m ,double* A, bool ascending)
{

	int n=m/2;

	double a1, a2;

	for(int i=0;i<n;i++)
	{

		if(A[i]>A[i+n])
		{

			a1=A[i+n];
			a2=A[i];

		} else {

			a1=A[i];
			a2=A[i+n];

		}

		A[i]   = ( ascending ? a1 : a2 );
		A[i+n] = ( ascending ? a2 : a1 );

	}
	return 0;

}

int bitonic_sort_bitonic(int n, double* C, bool ascending)
{

	int k=1;
	while(n>1)
	{

		for(int i=0;i<k;i++)
		{

			bitonic_split_sequential(n,&(C[i*n]),ascending);

		}
		k=k*2;

		n=n/2;

	}
	return 0;

}

int bitonic_sort_sequential_old(int n, double* A)
{

	int k=1;

	while(n>1)
	{
		k=k*2;
		n=n/2;
		for(int i=0;i<n;i++)
		{

			bitonic_split_sequential(k,&(A[i*k]),(i+1) % 2);
		}
		for(int j=0;j<8;j++)
		{
			cout << A[j] << " ";
		}
		cout << endl;



	}
	return 0;

}

int bitonic_sort_sequential(int n, double* A, bool ascending)
{

	

	return 0;	

}


int main () {
	double A[8]={17.0, 5.0, 9.0, 6.0, 11.0, 1.0, 3.0, 7.0};
	bitonic_sort_sequential(8,A);

	for(int j=0;j<8;j++)
	{
		//		cout << A[j] << " ";
	}
	//	cout << endl;

}
