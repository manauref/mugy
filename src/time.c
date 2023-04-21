/* mugy: time.c
 *
 * Functions to manipulate the time struct and do time integration.
 *
 */
#include "mh_time.h"
#include "mh_io_utilities.h"
#include "mh_dydt.h"
#include <stdlib.h>  // for malloc.
#include <stdio.h>  // for sprintf.
#include <math.h>  // for round.

struct mugy_time *mugy_time_alloc() {
  // Allocate the time object.
  struct mugy_time *time = (struct mugy_time *) malloc(sizeof(struct mugy_time));
  return time;
}

void mugy_time_wcstamp(struct timespec *wcstamp) {
  // Record the wall-clock time in a wc stamp.
  clock_gettime(CLOCK_MONOTONIC, wcstamp);
};

double mugy_time_sec(struct timespec wcstamp) {
  // Return the time of a wall-clock time in seconds.
  return wcstamp.tv_sec + 1e-9*wcstamp.tv_nsec;
}

struct timespec mugy_time_diff(struct timespec start, struct timespec end) {
  // Compute the time difference between two wall-clock stamps.
  struct timespec tm;
  if ((end.tv_nsec - start.tv_nsec)<0) {
    tm.tv_sec  = end.tv_sec-start.tv_sec-1;
    tm.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    tm.tv_sec  = end.tv_sec-start.tv_sec;
    tm.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return tm;
}

double mugy_time_elapsed_sec(struct timespec prev) {
  // Return the time elapsed since the a previous wall-clock stamp (prev).
  struct timespec curr;
  mugy_time_wcstamp(&curr);
  return mugy_time_sec(mugy_time_diff(prev, curr));
}

char *mugy_time_datetime() {
  // Obtain a string with the current date and time.
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  static char out[17];
  sprintf(out, "%d-%02d-%02d %02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  return out;
}

void mugy_time_init(struct mugy_time *time, mint world_rank, bool isRestart) {
  // Initialize time stepping.
  time->dt      = time->pars.dt;
  time->dt_prev = time->dt;
  time->dt_init = time->dt;
  time->dt_next = time->dt;
  time->dt_max  = time->dt*100.;  // Maximum allowed time step.
  time->ttol    = time->dt*1.e-2;

  time->simTime   = 0.;
  time->framesOut = 0;

  valPrintS_real(time->dt, "\n Initial time step: dt =", "\n", world_rank);

  // Time rate at which to output, adjust hyperdiffusion and adjust time step.
  time->tRateOutput = time->pars.endTime/((double) time->pars.nFrames);

  // Set the time rate at which to enter messages in the log and
  // determine the number of entries in log files/screen so far.
  time->tRateLogEntry = 1.e-2*time->pars.endTime;
  if (isRestart)
    time->logEntries = round(time->simTime/time->tRateLogEntry);
  else
    time->logEntries = 0;

}

void mugy_time_step_euler(mint outIdx, mint inIdx, mint dotIdx,
  enum mugy_op_types op, double time, double dt, struct mugy_population *pop) {
  // Take an Euler step of size dt.
  // Assume that:
  //   - The time-rate-of-change is in the step index momk[dotIdx].
  //   - The current state is in step index momk[inIdx].
  //   - The new state will be placed in in step index momk[outIdx].
  struct mugy_array *momIn  = pop->local->momk[inIdx];
  struct mugy_array *momOut = pop->local->momk[outIdx];
  struct mugy_array *momDot = pop->local->momk[dotIdx];

#ifdef USE_GPU
  enum mugy_resource_calc onResource = MUGY_DEVICE_CALC;
#else
  enum mugy_resource_calc onResource = MUGY_HOST_CALC;
#endif

  if (op == MUGY_OP_ASSIGN)
    mugy_array_axpy_assign(momOut, dt, momDot, momIn, onResource);
  else if (op == MUGY_OP_INCREMENT)
    mugy_array_axpy_increment(momOut, dt, momDot, momIn, onResource);
}

void mugy_time_step_rk4(double time, double dt, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft) {
  // Step the solution forward by one time step of size dt
  // using 4th-order Runge Kutta.
  
#ifdef USE_GPU
  enum mugy_resource_calc onResource = MUGY_DEVICE_CALC;
#else
  enum mugy_resource_calc onResource = MUGY_HOST_CALC;
#endif

  // The following is a convoluted way of doing RK4 so that only 3 storage locations.
  // MF 2023/04/19: this may not be the official way of doing low-storage RK4,
  //                it is simply something I devised by looking at the Butcher tableau.
  double tnew = time;
  mugy_dydt(1, 0, tnew, pop, field, grid, fft);
  mugy_time_step_euler(1, 0, 0.5*dt, 1, MUGY_OP_ASSIGN   , tnew, pop);
  mugy_field_poisson_solve(field, pop, grid, 1);

  tnew = time+0.5*dt;
  mugy_dydt(2, 1, tnew, pop, field, grid, fft);
  mugy_time_step_euler(1, 0,     dt, 2, MUGY_OP_INCREMENT, tnew, pop);
  mugy_time_step_euler(2, 0, 0.5*dt, 2, MUGY_OP_ASSIGN   , tnew, pop);
  mugy_field_poisson_solve(field, pop, grid, 2);

  tnew = time+0.5*dt;
  mugy_dydt(2, 2, tnew, pop, field, grid, fft);
  mugy_time_step_euler(2, 0,     dt, 2, MUGY_OP_ASSIGN   , tnew, pop);
  mugy_field_poisson_solve(field, pop, grid, 2);
  mugy_array_increment(pop->local->momk[1], pop->local->momk[2], onResource);

  tnew = time+dt;
  mugy_dydt(2, 2, tnew, pop, field, grid, fft);
  mugy_time_step_euler(1, 1, 0.5*dt, 2, MUGY_OP_ASSIGN   , tnew, pop);

  mugy_array_ax_assign(pop->local->momk[0], 1./3., pop->local->momk[1], onResource);
  mugy_field_poisson_solve(field, pop, grid, 0);

//  double tnew = time;
//  mugy_dydt(1, 0, tnew, pop, field, grid, fft);
//  mugy_time_step_euler(1, 0, dt, 1, MUGY_OP_ASSIGN   , tnew, pop);
//  mugy_field_poisson_solve(field, pop, grid, 1);  // Solve Poisson to get phik.
//  mugy_array_ax_assign(pop->local->momk[0], 1., pop->local->momk[1], onResource);
};

void mugy_time_advance(double time, double dt, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft) {
  // Step the solution forward by one time step of size dt.
  
#if TIME_STEPPER==4
  // Step with 4th-order Runget Kutta.
  mugy_time_step_rk4(time, dt, pop, field, grid, fft);
#endif
};

void mugy_time_free(struct mugy_time *time) {
  // Free memory associated with time object.
  free(time);
}
