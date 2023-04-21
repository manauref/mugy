/* mugy: mh_dydt_dev.h
 *
 * Compute the time rate of change of the moment vector on the device.
 *
 */
#pragma once

#include "mh_dydt.h"

// Compute <J_0>phik, 0.5*hatLap<J_0>phik, 1+0.5*hatLap+hathatLap)<J_0>phik on device.
void mugy_dydt_phikFLR_dev(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid);

// Apply the linear terms on the device.
void mugy_dydt_linear_dev(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid);

// Compute the time rate of change dy/dt, where y is the vector of moments.
void mugy_dydt_dev(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft);
