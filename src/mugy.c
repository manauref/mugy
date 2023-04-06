/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  struct mugy_grid gridG, gridL;
  struct mugy_timeSetup timePars;
  struct mugy_population popG, popL;
  struct mugy_fieldParameters fieldPars; 
  struct mugy_ioManager ioMan;
  struct mugy_ffts fftMan;
  struct mugy_timeState tState;

  init_mpi(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n" );

  // Run the full initialization.
  init_all(argc, argv, &ioMan, &gridG, &gridL, &timePars, &popG, &popL, &fieldPars, &fftMan);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  fft_terminate(&fftMan);  // Deallocate FFT memory.
  io_terminate(&ioMan);  // Close IO interface. Need to do it before freeing fields.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  terminate_mpi();  // Finalize MPI.
//  terminate_all();  // Call termination of all parts of mugy.

  return 0;
}
