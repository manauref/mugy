/* mugy: mh_timestepping.h
 *
 * Structs and functions in time integration.
 *
 */
#ifndef MUGY_TIMESTEPPING
#define MUGY_TIMESTEPPING

#include "mh_data.h"

struct mugy_timeSetup {
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

struct mugy_timeState {
  real simTime;
  mint time;
  mint framesOut;
  mint hdAdjusts;
  mint dtAdjusts;
  real dt;
};

#endif
