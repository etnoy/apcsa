#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


inline float rastrigin(double* x, int n) {
  float sum = 0;
  int i;
  for(i = 0; i < n; i++) {
    sum += x[i] * x[i] - 10 * cos(2.0*M_PI*x[i]);
  }
  return 10*n+sum;
}

void copyState(double* from, double* to, int n) {
  int i;
  for(i = 0; i < n; i++) {
    to[i] = from[i];
  }
}

void printState(double* aState, int n) {
  int i;
  for(i = 0; i < n; i++) {
    printf("x[%d]=%2.10f\n", i, aState[i]);
  }
}

int main(int argc, char** argv) {
  ///////////////////////////////////////////////////////////////////
  // Initialize the parameters
  ///////////////////////////////////////////////////////////////////
  unsigned int dim = 3;
  unsigned int lambda = 20; //assume is a multiple of the team-size
  unsigned int max_generation = 100;
  unsigned int number_of_groups = 5; // of course smaller than num procs
  double variance;  
  int bbox[dim];
  
  int i;
  for(i = 0; i < dim; i++) {
    bbox[i] = 5;
  }

  // Compute the machine epsilon
  float epsilon = 1.0f;
  do {
    epsilon /= 2.0f;
  }
  while ((float)(1.0 + (epsilon/2.0f)) != 1.0f);
  epsilon *= 4;
  
  ///////////////////////////////////////////////////////////////////
  //initialize mpi variables
  ///////////////////////////////////////////////////////////////////
  int rank, size, ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);

  ///////////////////////////////////////////////////////////////////
  // TODO: 
  // - determine the size of the groups
  // - build the group communicators
  // - assign team_rank to each process (the rank within the group)
  ///////////////////////////////////////////////////////////////////
  int size_of_groups = ;
  int color = ; //the same as 'team'
  MPI_Comm team_comm;
  int team_rank, team_size;
  

  ///////////////////////////////////////////////////////////////////
  // Initialization of the random number generator
  // TODO:
  // Set the seed of the random number generator with 
  // gsl_rng_set(const gsl_rng * r, unsigned long int s)
  // The passionate may use the system time
  ///////////////////////////////////////////////////////////////////
  gsl_rng* rand_gen = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(...);

  ///////////////////////////////////////////////////////////////////
  // TODO: For each team, determine the value of its variance.
  ///////////////////////////////////////////////////////////////////
  variance = ...;
  
  ///////////////////////////////////////////////////////////////////
  // The init of state vectors
  ///////////////////////////////////////////////////////////////////
  double x_parent[dim]; //equals x0
  for(i = 0; i < dim; i++){
    x_parent[i] = (gsl_rng_uniform(rand_gen) - 0.5) * 2.0 * bbox[i];
  }
  float f_x_parent = rastrigin(x_parent, dim);

   // initialize the global best and the best of generation x
  double x_best_generation[dim];
  float f_x_best_generation = FLT_MAX;
  copyState(x_parent,x_best_generation,dim);

  // every process stores the overall best of the team
  double x_best_team[dim];
  float f_x_best_team = FLT_MAX;

  ///////////////////////////////////////////////////////////////////
  // start the evolution
  ///////////////////////////////////////////////////////////////////
  unsigned int generationIt;
  for(generationIt = 0; generationIt < max_generation; generationIt++) {
    f_x_best_generation = FLT_MAX;
    unsigned int success = 0;
    double x[dim];
    int lambdaIt;

  ///////////////////////////////////////////////////////////////////
  // start sampling the offspring
  ///////////////////////////////////////////////////////////////////
    for(lambdaIt = 0; lambdaIt < lambda/team_size; lambdaIt++) {

      // generate the new sample with the actual variance
      for(i = 0; i < dim; i++){
        x[i] = x_parent[i] + gsl_ran_gaussian(rand_gen,variance);
      }

      float f_x = rastrigin(x, dim);

      // check if the sample is better than the current best sample.
      // if yes, copy the values to the generations best...
      if(f_x < f_x_best_generation) {
        f_x_best_generation = f_x;
        copyState(x, x_best_generation, dim);
      }
    }

    ///////////////////////////////////////////////////////////////////
    // TODO:
    // figure out which process found the best value:
    // 1) The root of each team gathers all f(x) the processes found
    // 2) Root finds the index of the best process and broadcasts this within the group.
    // 3) This process then broadscasts its result (with the state vec)
    ///////////////////////////////////////////////////////////////////
    float best_team_value;
    float* all_best_values = ;    
    int best_rank_i;




    // after having x_bes_generation of all the team members, select 
    // the best offspring as the new parent.
    copyState(x_best_generation, x_parent, dim);
    f_x_parent = f_x_best_generation;

    // store the overall best result of this team:
    if(f_x_best_generation < f_x_best_team) {
      f_x_best_team = f_x_best_generation;
      copyState(x_best_generation, x_best_team, dim);
    }
    // check if we already found the global minimum - will be seldom
    if(f_x_parent < epsilon) {
      break;
    }

    free(all_best_values);
  }


  ///////////////////////////////////////////////////////////////////
  // TODO: Collect the results:
  // - All the roots of the teams send the best team result (f(x)), 
  //   x, variance and team number to root.
  // - root searches for the best result and prints it to the screen.
  ///////////////////////////////////////////////////////////////////
  
 
  
  // free memory
  gsl_rng_free(rand_gen);

  // finalize MPI
  ierr = MPI_Finalize();
  return 0;
}
