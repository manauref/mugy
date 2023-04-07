/*
   mugy
   A reduced gyrofluid code for ITG-ETG multiscale simulations.
*/

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  struct mugy_grid gridG, gridL;
  struct mugy_comms comms;
  struct mugy_timeSetup timePars;
  struct mugy_population popG, popL;
  struct mugy_fieldParameters fieldPars; 
  struct mugy_ioManager ioMan;
  struct mugy_ffts fftMan;
  struct mugy_timeState tState;

  comms_init(&comms, argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n", comms.world.rank );

  // Run the full initialization.
  init_all(argc, argv, &comms, &ioMan, &gridG, &gridL, &timePars, &popG, &popL, &fieldPars, &fftMan);

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  io_terminate(&ioMan);  // Close IO interface. Need to do it before freeing fields.
  fft_terminate(&fftMan);  // Deallocate FFT memory.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  comms_terminate(&comms);  // Finalize communications.
  terminate_all();  // Call termination of all parts of mugy.

  return 0;
}
