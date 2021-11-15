/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mugy.h"

int main(int argc, char *argv[]) {

  struct grid gridG, gridL;
  struct timeSetup timePars;
  struct speciesParameters specParsG, specParsL;
  struct ioSetup myIO;
  struct fieldParameters fieldPars; 

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &myIO, &gridG, &timePars, &specParsG, &fieldPars);
  init_comms(gridG, specParsG);

  // Set the number of cells in Fourier space and aliased real space.
  init_global_grids(&gridG);
  distributeDOFs(gridG, specParsG, &gridL, &specParsL);

  allocate_fields(gridL, specParsL);

  set_initialCondition(gridL, specParsL);

//  printf(" Number of time steps and frames:           Nt       =%8d   |  nFrames  =%6d\n", 1000, timePars.nFrames);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  free_fields();
  free_grid(&gridL);
  free_speciesPars(&specParsL);
  free_grid(&gridG);
  free_speciesPars(&specParsG);

  terminate_mpi();  // Finalize MPI.

  return 0;
}
