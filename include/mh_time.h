/* mugy: mh_timestepping.h
 *
 * Structs and functions in time integration.
 *
 */
#pragma once

#include "mh_data.h"
#include <time.h>

// Time-related parameters read from input file.
struct mugy_time_pars {
  real dt;               // Time step.
  real endTime;          // Absolute end time (from t=0 of first simulation).
  mint nFrames;          // Absolute frames to output (from t=0 of first simulation).
  mint ark_kySplit;      // multirate splitting index: highest "slow" wavenumber (higher are "fast").
  mint ark_fastTableExp; // Butcher table index for fast explicit method.
  mint ark_fastTableImp; // Butcher table index for fast implicit method.
  mint ark_slowTable;    // Butcher table index for slow method.
  real ark_dtFast;       // fixed 'fast' time step to use.
  real ark_rtol;         // relative solution tolerance for temporal adaptivity.
  real ark_atol;         // absolute solution tolerance for temporal adaptivity
  mint ark_ewtScaling;   // which error weight scaling to use.
};

// Wall-clock stamps, for timing various components.
struct mugy_time_wcstamps {
  struct timespec all;       // Whole program.
  struct timespec init;      // Initialization.
  struct timespec timeloop;  // Main time loop.
  struct timespec fin;       // Finalization.
};

struct mugy_time {
  struct mugy_time_pars pars;
  struct mugy_time_wcstamps wcs;
  double simTime;    // Simulation time.
  double dt;         // Time step size.
  long steps;      // Number of steps.
  long framesOut;  // Number of frames.
  long hdAdjusts;
  long dtAdjusts;
  // We need several other time steps; previous, initial, next and max.
  double dt_prev, dt_init, dt_next, dt_max;
  double ttol;  // Tolerance used to trigger actions in time loop.
  double tRateOutput;  // Time rate at which to output.
};

struct mugy_time *mugy_time_alloc();

// Record the wall-clock time in a wc stamp.
void mugy_time_wcstamp(struct timespec *wcstamp);
// Return the time of a wall-clock time in seconds.
double mugy_time_sec(struct timespec wcstamp);
// Compute the time difference between two wall-clock stamps.
struct timespec mugy_time_diff(struct timespec start, struct timespec end);
// Return the time elapsed since the a previous wall-clock stamp (prev).
double mugy_time_elapsed_sec(struct timespec prev);
// Obtain a string with the current date and time.
char *mugy_time_datetime();

// Initialize time stepping.
void mugy_time_init(struct mugy_time *time, mint world_rank);

// Free memory associated with time object.
void mugy_time_free(struct mugy_time *time);
