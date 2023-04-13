/* mugy
 *
 * A reduced gyrofluid code for ITG-ETG multiscale simulations.
 *
 */

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  struct mugy_timeSetup timePars;
  struct mugy_fieldParameters fieldPars; 
  
  struct mugy_timeState tState;

  // ............ INITIALIZATION ............ //

  struct mugy_comms *comms = mugy_comms_init(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n", comms->world.rank );

  struct mugy_ioManager *ioMan = mugy_io_init(*comms);  // Initialize IO interface.

  struct mugy_grid *grid = mugy_grid_alloc();  // Allocate grid object.
  struct mugy_population *pop = mugy_population_alloc();  // Allocate population object.

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &ioMan->setup, grid, &timePars, pop, &fieldPars, comms->world.rank);

#ifdef USE_GPU
  // Initialize devices (GPUs) if any.
  device_init(comms);
#endif

  mugy_comms_sub_init(comms, *grid, *pop);  // Initialize sub-communicators.

  // Set the number of cells in Fourier space and aliased real space.
  mugy_grid_init_global(grid, comms->world.rank);

  // Decompose the x,y,z,s domains amongst MPI processes.
  mugy_comms_distributeDOFs(*comms, grid, pop);

  mugy_population_allocate_moments(pop, *grid);  // Allocate dynamic moments.

  struct mugy_ffts *fftMan = mugy_fft_init(*grid, pop->local, *comms);  // Initialize FFT infrastructure.

  mugy_io_setup_files(ioMan, *grid, *pop);  // Setup IO files.

  set_initialCondition(*grid, pop, fftMan, ioMan);  // Impose ICs.

  // Write initial conditions.
  mugy_io_write_mugy_array(ioMan, "momk", NULL, *pop->local.momk);

  MPI_Barrier(comms->world.comm); // To avoid starting time loop prematurely.
  // ............ END OF INITIALIZATION ............ //



  // ............ DEALLOCATE / FINALIZE ............ //
  MPI_Barrier(comms->world.comm); // To avoid premature deallocations.
  mugy_io_terminate(ioMan);  // Close IO interface. Need to do it before freeing fields.
  mugy_fft_terminate(fftMan);  // Deallocate FFT memory.

  mugy_grid_free(grid);
  mugy_population_free(pop);

  mugy_comms_terminate(comms);  // Finalize communications.

  // ............ END OF DEALLOCATE / FINALIZE ............ //
  return 0;
}
