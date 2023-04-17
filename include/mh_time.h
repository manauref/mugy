/* mugy: mh_timestepping.h
 *
 * Structs and functions in time integration.
 *
 */
#pragma once

#include "mh_data.h"

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

struct mugy_time {
  struct mugy_time_pars pars;
  real time;       // Simulation time.
  mint framesOut;
  mint hdAdjusts;
  mint dtAdjusts;
  real dt;         // Time step size.
};

struct mugy_time *mugy_time_alloc();

void mugy_time_free(struct mugy_time *time);
