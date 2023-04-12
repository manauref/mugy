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
  
  struct mugy_timeState tState;

  // ............ INITIALIZATION ............ //

  struct mugy_comms *comms = comms_init(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n", comms->world.rank );

  struct mugy_ioManager *ioMan = io_init(*comms);  // Initialize IO interface.

    // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &ioMan->setup, &gridG, &timePars, &popG, &fieldPars, comms->world.rank);

#ifdef USE_GPU
  // Initialize devices (GPUs) if any.
  device_init(comms);
#endif

  comms_sub_init(comms, gridG, popG);  // Initialize sub-communicators.

  // Set the number of cells in Fourier space and aliased real space.
  init_global_grids(&gridG, comms->world.rank);

  // Decompose the x,y,z,s domains amongst MPI processes.
  distributeDOFs(*comms, gridG, popG, &gridL, &popL);

  allocate_dynfields(gridL, &popL);  // Allocate dynamic fields.

  struct mugy_ffts *fftMan = mugy_fft_init(gridG, gridL, popL, *comms);  // Initialize FFT infrastructure.

  setup_files(ioMan, gridG, gridL, popG, popL);  // Setup IO files.

  set_initialCondition(gridG, gridL, popG, &popL, fftMan, ioMan);  // Impose ICs.

  // ............ END OF INITIALIZATION ............ //

  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  io_terminate(ioMan);  // Close IO interface. Need to do it before freeing fields.
  mugy_fft_terminate(fftMan);  // Deallocate FFT memory.

  free_fields();
  free_grid(&gridL);
  free_population(&popL);
  free_grid(&gridG);
  free_population(&popG);

  comms_terminate(comms);  // Finalize communications.
  terminate_all();  // Call termination of all parts of mugy.

  return 0;
}
