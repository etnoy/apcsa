#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<string.h>

unsigned int MATRIX_SIZE = 100;

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
  double* my_A; // this is calculated in the beginning and remains const
  double* A;
  double* B;
  double* C;
} COMM_INFO;

void generate_communicators(COMM_INFO* comm_info) {
  int origrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &(origrank));
  MPI_Comm_size(MPI_COMM_WORLD, &(comm_info->size));

  unsigned int i;
  int q = (int)sqrt(comm_info->size);
  comm_info->q = q;
  
  int dims[2] = {q, q};
  int periods[2] = {1, 1};
  MPI_Comm grid_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

  MPI_Comm_rank(grid_comm, &(comm_info->rank));

  int coords[2];
  MPI_Cart_coords(grid_comm, comm_info->rank, 2, coords);
  comm_info->row_index = coords[0];
  comm_info->col_index = coords[1];
  comm_info->row_start = MATRIX_SIZE/q * comm_info->row_index;
  comm_info->col_start = MATRIX_SIZE/q * comm_info->col_index;

  int free_coords[2];
  free_coords[0] = 0;
  free_coords[1] = 1;
  MPI_Cart_sub(grid_comm, free_coords, &(comm_info->row_comm));
  free_coords[0] = 1;
  free_coords[1] = 0;
  MPI_Cart_sub(grid_comm, free_coords, &(comm_info->col_comm));

  MPI_Comm_rank(comm_info->col_comm, &(comm_info->col_comm_rank));
  MPI_Comm_rank(comm_info->row_comm, &(comm_info->row_comm_rank));
  
}

double* generate_sub_matrix(COMM_INFO* comm_info) {
  double* sub_matrix;
  int sub_matrix_size = MATRIX_SIZE / comm_info->q;

  sub_matrix = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));
  int j,i;
  for(i = 0; i < sub_matrix_size; i++) {
    for(j = 0; j < sub_matrix_size; j++) {
      sub_matrix[i*sub_matrix_size + j] =
        cos((comm_info->row_start + i) * MATRIX_SIZE
            + comm_info->col_start + j);
    }
  }
  return sub_matrix;
}
void collect_and_print_element(COMM_INFO* comm_info,
                               unsigned int i, unsigned int j) {
  int sub_matrix_size = MATRIX_SIZE / comm_info->q;
  double element_ij;

  if(i < sub_matrix_size && j < sub_matrix_size) {//no communication needed!
    element_ij = comm_info->C[i * sub_matrix_size + j];
  } else {
    if(comm_info->rank == 0) {
      MPI_Status stat;
      MPI_Recv(&element_ij, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
               9, MPI_COMM_WORLD, &stat);
    }

    else if((comm_info->row_start <= i &&
             i < comm_info->row_start + sub_matrix_size &&
             comm_info->col_start <= j &&
             j < comm_info->col_start + sub_matrix_size)) {
      int my_i = i - comm_info->row_start;
      int my_j = j - comm_info->col_start;
      MPI_Send(&(comm_info->C[my_i*sub_matrix_size + my_j]),
               1, MPI_DOUBLE,
               0, 9, MPI_COMM_WORLD);
    }
  }

  if(comm_info->rank == 0) {
    printf("Element C(%d,%d) = %lf.\n", i, j, element_ij);
  }
}

int multiply_local(COMM_INFO* comm_info) {
  int i,j, r;
  int s = MATRIX_SIZE / comm_info->q;

  double* A = comm_info->A;
  double* B = comm_info->B;
  double* C = comm_info->C;

  for(i = 0; i < s; i++) {
    for(j = 0; j < s; j++) {
      for(r = 0; r < s; r++) {
        C[i*s + j] +=
          (A[i * s + r]) * (B[r * s + j]);
      }
    }
  }
}

int main(int argc, char** argv) {

  COMM_INFO comm_info;
  MPI_Init(&argc, &argv);
  generate_communicators(&comm_info);
  int q = comm_info.q;
  int sub_matrix_size = MATRIX_SIZE / q;

  comm_info.my_A = generate_sub_matrix(&comm_info);
  comm_info.B = generate_sub_matrix(&comm_info);
  comm_info.C = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));
  comm_info.A = malloc(sub_matrix_size * sub_matrix_size * sizeof(double));

  memset(comm_info.C, 0.0, sub_matrix_size*sub_matrix_size * sizeof(double));

  ///////////////////////////////////////////////////////////////////
  // TODO: subquestion e)
  ///////////////////////////////////////////////////////////////////

  int rank_to_send_to = (comm_info.col_comm_rank + q - 1) % q;
  int rank_to_get_from = (comm_info.col_comm_rank + 1) % q;

  MPI_Status stat;
  int k;

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for(k = 0; k < q; k++) {
    int k_bar = (k + comm_info.row_index) % q;
    //separate sender and recv, we don't want to overwrite my_A;
    printf("%d %d\n", k_bar, comm_info.row_comm_rank);
    if(k_bar == comm_info.row_comm_rank) {

      fflush(stdout);
      MPI_Bcast(comm_info.my_A,
                sub_matrix_size * sub_matrix_size,
                MPI_DOUBLE,
                k_bar,
                comm_info.row_comm);

      //!!!! take care, also A has to have the correct value!
      memcpy(comm_info.A, comm_info.my_A,
             sub_matrix_size * sub_matrix_size * sizeof(double));

    } else {

      fflush(stdout);
      MPI_Bcast(comm_info.A,
                sub_matrix_size * sub_matrix_size,
                MPI_DOUBLE,
                k_bar,
                comm_info.row_comm);

    }

    multiply_local(&comm_info);
    MPI_Sendrecv_replace(comm_info.B,
                         sub_matrix_size * sub_matrix_size,
                         MPI_DOUBLE,
                         rank_to_send_to,
                         comm_info.col_comm_rank+q*k,
                         rank_to_get_from,
                         rank_to_get_from+q*k,
                         comm_info.col_comm,
                         &stat);
  }
  collect_and_print_element(&comm_info,98,99);

  free(comm_info.A);
  free(comm_info.B);
  free(comm_info.C);
  free(comm_info.my_A);
  MPI_Finalize();

}
