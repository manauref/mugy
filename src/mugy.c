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
  struct timeState tState;

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Run the full initialization.
  init_all(argc, argv, &myIO, &gridG, &gridL, &timePars, &popG, &popL, &fieldPars);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  terminate_io();  // Close IO interface. Need to do it before freeing fields.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  terminate_all();  // Call termination of all parts of mugy.

  return 0;
}
