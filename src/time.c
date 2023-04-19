/* mugy: time.c
 *
 * Functions to manipulate the time struct and do time integration.
 *
 */
#include "mh_time.h"
#include "mh_io_tools.h"
#include <stdlib.h>  // for malloc.
#include <stdio.h>  // for sprintf.

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

void mugy_time_init(struct mugy_time *time, mint world_rank) {
  // Initialize time stepping.
  time->dt      = time->pars.dt;
  time->dt_prev = time->dt;
  time->dt_init = time->dt;
  time->dt_next = time->dt;
  time->dt_max  = time->dt*100.;  // Maximum allowed time step.
  time->ttol    = time->dt*1.e-2;

  time->simTime   = 0.;
  time->framesOut = 0;

  // Time rate at which to output, adjust hyperdiffusion and adjust time step.
  time->tRateOutput = time->pars.endTime/((double) time->pars.nFrames);

  valPrintS_real(time->dt, "\n Initial time step: dt =", "\n", world_rank);
}

void mugy_time_free(struct mugy_time *time) {
  // Free memory associated with time object.
  free(time);
}
