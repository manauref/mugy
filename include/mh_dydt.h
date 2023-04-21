/* mugy: mh_dydt.h
 *
 * Compute the time rate of change of the moment vector.
 *
 */
#pragma once

#include "mh_macros.h"
#include "mh_population.h"
#include "mh_field.h"
#include "mh_grid.h"
#include "mh_ffts.h"

// Compute <J_0>phik, 0.5*hatLap<J_0>phik, 1+0.5*hatLap+hathatLap)<J_0>phik.
void mugy_dydt_phikFLR(struct mugy_population *pop, struct mugy_field *field, struct mugy_grid *grid);

// Apply the linear terms.
void mugy_dydt_linear(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid);

// Compute the time rate of change dy/dt, where y is the vector of moments.
void mugy_dydt(mint outIdx, mint inIdx, double time, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_grid *grid, struct mugy_fft *fft);
