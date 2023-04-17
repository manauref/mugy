/* mugy: mh_flr.h
 *
 * Create various FLR (finite Larmor radius) operators.
 *
 */
#pragma once

#include "mh_population.h"
#include "mh_grid.h"

struct mugy_flr {
  struct mugy_array *Gamma0   ;  // Exponentially weighted 0th modified Bessel function of the 1st kind.
  struct mugy_array *avgJ0    ;  // <J_0>=Gamma0^{1/2}.
  struct mugy_array *hatLap   ;  // \hat{nabla}_perp.
  struct mugy_array *hathatLap;  // \hat{hat{nabla}}_perp.
};

struct mugy_flr *mugy_flr_init(struct mugy_population *pop, struct mugy_grid *grid);

void mugy_flr_free(struct mugy_flr *flr);
