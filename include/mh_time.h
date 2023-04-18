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
  real time;       // Simulation time.
  mint framesOut;
  mint hdAdjusts;
  mint dtAdjusts;
  real dt;         // Time step size.
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

// Free memory associated with time object.
void mugy_time_free(struct mugy_time *time);
