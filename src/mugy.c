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

  mugy_time_wcstamp(&time->wcs.init);  // Start timing initialization.

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

  // Initialize time independent operators.
  mugy_constop_init(pop, grid, field, io);

  // Impose ICs.
  set_initialConditions(pop, field, grid, fft, io);

  // Write initial conditions.
  mugy_io_write_mugy_array(io, "momk", NULL, pop->local->momk[0]);

  // Initialize time stepping infrastructure.
  mugy_time_init(time, comms->world->rank);

  MPI_Barrier(comms->world->comm); // Avoid starting time loop prematurely.

  double tm_init = mugy_time_elapsed_sec(time->wcs.init);  // Finish timing initialization.
  valPrintS_real(tm_init, "\n Time spent on initialization: ", " s\n", comms->world->rank); 
  // ............ END OF INITIALIZATION ............ //


  // ............ TIME LOOP ............ //
  mugy_time_wcstamp(&time->wcs.timeloop);  // Start timing time loop.
  valPrintS_char(mugy_time_datetime(), "\n Entering time loop on ", "\n", comms->world->rank); 

  bool continue_stepping = true;
  while (continue_stepping) {

    time->dt = time->dt_next;

    time->steps   += 1;         // Time steps taken.
    time->simTime += time->dt;  // Current simulation time.

    double tOutNext = ((double) (time->framesOut+1))*time->tRateOutput;
    if ( (fabs(time->simTime - tOutNext) <= time->ttol) ||
         (fabs(time->simTime - tOutNext) < fabs(time->simTime+time->dt_next - tOutNext)) ) {
//      isNaNorInf(pop, field, grid);    // Check if solutions diverged.

      time->framesOut = time->framesOut+1;

//      mugy_io_write_dynamic_diagnostics(pop, field, grid);    // Append evolving quantities to output files.
//
//      mugy_io_write_restart();     // Write file used for restarts.
    }

    if (fabs(time->simTime-time->pars.endTime) <= time->ttol) {
      // Time loop completed. Punch out.
      continue_stepping = false;
      break;
    } else if (time->simTime+time->dt_next > time->pars.endTime) {
      // The end is less than a time step away. Take a smaller time step.
      time->dt_next = time->pars.endTime - time->simTime;
    }
  }

  double tm_tloop = mugy_time_elapsed_sec(time->wcs.timeloop);  // Finish timing time loop.
  valPrintS_real(tm_tloop, "\n Time spent on time loop: ", " s\n", comms->world->rank); 
  // ............ END OF TIME LOOP ............ //


  // ............ DEALLOCATE / FINALIZE ............ //
  mugy_time_wcstamp(&time->wcs.fin);  // Start timing termination.
  MPI_Barrier(comms->world->comm); // Avoid premature deallocations.
  mugy_io_terminate(io);  // Close IO interface BEFORE freeing arrays written in time loop.
  mugy_fft_terminate(fft);
  mugy_time_free(time);
  mugy_field_free(field);
  mugy_grid_free(grid);
  mugy_population_free(pop);

  // Need to put this before subcommunicators get deallocated.
  double tm_fin = mugy_time_elapsed_sec(time->wcs.fin);  // Finish timing termination.
  valPrintS_real(tm_fin, "\n Time spent on termination: ", " s\n", comms->world->rank); 

  mugy_comms_terminate(comms);  // Finalize communications.
  // ............ END OF DEALLOCATE / FINALIZE ............ //

  return 0;
}
