/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mugy.h"

int main(int argc, char *argv[]) {

  struct grid gridG, gridL;
  struct timeSetup timePars;
  struct population popG, popL;
  struct ioSetup myIO;
  struct fieldParameters fieldPars; 

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &myIO, &gridG, &timePars, &popG, &fieldPars);
  init_comms(gridG, popG);

  // Set the number of cells in Fourier space and aliased real space.
  init_global_grids(&gridG);
  distributeDOFs(gridG, popG, &gridL, &popL);

  allocate_fields(gridL, popL);

  set_initialCondition(gridL, popL);

//  printf(" Number of time steps and frames:           Nt       =%8d   |  nFrames  =%6d\n", 1000, timePars.nFrames);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  terminate_mpi();  // Finalize MPI.

  return 0;
}
