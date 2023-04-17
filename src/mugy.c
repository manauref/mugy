/* mugy
 *
 * A reduced gyrofluid code for ITG-ETG multiscale simulations.
 *
 */

#include "mh_mugy.h"

int main(int argc, char *argv[]) {

  // ............ INITIALIZATION ............ //
  struct mugy_comms *comms = mugy_comms_init(argc, argv);  // Initialize MPI interface.

  r0printf("\n     --> Welcome to mugy <--    \n\n", comms->world->rank );

  // Initialize IO interface.
  struct mugy_io *io = mugy_io_init(comms);

  // Allocate grid, population, field, time objects.
  struct mugy_grid *grid      = mugy_grid_alloc();
  struct mugy_population *pop = mugy_population_alloc();
  struct mugy_field *field    = mugy_field_alloc();
  struct mugy_time *time      = mugy_time_alloc();

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &io->setup, grid, &time->pars, pop, field, comms->world->rank);

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
  struct mugy_fft *fft = mugy_fft_init(grid, pop->local, comms);

  // Setup IO files.
  mugy_io_setup_files(io, grid, pop);

  // Initialize field object.
  mugy_field_init(field, grid, pop);

  // Initialize FLR and linear operators.
  mugy_flr_init(pop, grid, field, io);

  // Impose ICs.
  set_initialConditions(pop, field, grid, fft, io);

  // Write initial conditions.
  mugy_io_write_mugy_array(io, "momk", NULL, pop->local->momk[0]);

  MPI_Barrier(comms->world->comm); // Avoid starting time loop prematurely.
  // ............ END OF INITIALIZATION ............ //



  // ............ DEALLOCATE / FINALIZE ............ //
  MPI_Barrier(comms->world->comm); // Avoid premature deallocations.
  mugy_io_terminate(io);  // Close IO interface BEFORE freeing arrays written in time loop.
  mugy_fft_terminate(fft);
  mugy_time_free(time);
  mugy_field_free(field);
  mugy_grid_free(grid);
  mugy_population_free(pop);
  mugy_comms_terminate(comms);  // Finalize communications.
  // ............ END OF DEALLOCATE / FINALIZE ............ //
  return 0;
}
