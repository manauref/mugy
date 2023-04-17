/* mugy: mh_nonlinear.h
 *
 * Tools to compute the nonlinear terms (Poisson brackets).
 *
 */
#pragma once

#include "mh_population.h"
#include "mh_grid.h"
#include "mh_flr.h"

void mugy_nonlinear_constop_init(struct mugy_population *pop, struct mugy_grid *grid, struct mugy_flr *flr);
