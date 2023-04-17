/* mugy
 *
 * A reduced gyrofluid code for ITG-ETG multiscale simulations.
 *
 */

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  struct mugy_timeSetup timePars;
  struct mugy_timeState tState;

  // ............ INITIALIZATION ............ //
  struct mugy_comms *comms = mugy_comms_init(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n", comms->world->rank );

  // Initialize IO interface.
  struct mugy_ioManager *ioMan = mugy_io_init(comms);

  // Allocate grid, population and field objects.
  struct mugy_grid *grid      = mugy_grid_alloc();
  struct mugy_population *pop = mugy_population_alloc();
  struct mugy_field *field    = mugy_field_alloc();

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &ioMan->setup, grid, &timePars, pop, field, comms->world->rank);

#ifdef USE_GPU
  // Initialize devices (GPUs) if any.
  device_init(comms);
#endif

  // Initialize sub-communicators.
  mugy_comms_sub_init(comms, grid, pop);

  // Set the number of cells in Fourier space and aliased real space.
  mugy_grid_init_global(grid, comms->world->rank);

  // Decompose the x,y,z,s domains amongst MPI processes.
  mugy_comms_distributeDOFs(comms, grid, pop);

  // Allocate local arrays.
  mugy_population_alloc_local(pop->local, grid->local);

  // Initialize FFT infrastructure.
  struct mugy_ffts *fftMan = mugy_fft_init(grid, pop->local, comms);

  // Setup IO files.
  mugy_io_setup_files(ioMan, grid, pop);

  // Initialize field object.
  mugy_field_init(field, grid, pop);

  // Initialize FLR and linear operators.
  mugy_flr_init(pop, grid, field, ioMan);

  // Impose ICs.
  set_initialConditions(pop, field, grid, fftMan, ioMan);

  // Write initial conditions.
  mugy_io_write_mugy_array(ioMan, "momk", NULL, pop->local->momk[0]);

  MPI_Barrier(comms->world->comm); // Avoid starting time loop prematurely.
  // ............ END OF INITIALIZATION ............ //



  // ............ DEALLOCATE / FINALIZE ............ //
  MPI_Barrier(comms->world->comm); // Avoid premature deallocations.
  mugy_io_terminate(ioMan);  // Close IO interface BEFORE freeing arrays written in time loop.
  mugy_fft_terminate(fftMan);
  mugy_field_free(field);
  mugy_grid_free(grid);
  mugy_population_free(pop);
  mugy_comms_terminate(comms);  // Finalize communications.
  // ............ END OF DEALLOCATE / FINALIZE ............ //
  return 0;
}
