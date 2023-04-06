/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  struct grid gridG, gridL;
  struct timeSetup timePars;
  struct population popG, popL;
  struct mugy_ioManager ioMan;
  struct fieldParameters fieldPars; 
  struct timeState tState;

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Run the full initialization.
  init_all(argc, argv, &ioMan, &gridG, &gridL, &timePars, &popG, &popL, &fieldPars);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  terminate_io(&ioMan);  // Close IO interface. Need to do it before freeing fields.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  terminate_all();  // Call termination of all parts of mugy.

  return 0;
}
