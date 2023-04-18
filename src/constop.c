/* mugy: constop.c
 *
 * Precompute time independent operators.
 *
 */
#include "mh_constop.h"
#include "mh_flr.h"
#include "mh_linear.h"
#include "mh_nonlinear.h"
#include <math.h>

void mugy_constop_init(struct mugy_population *pop, struct mugy_grid *grid, struct mugy_field *field, struct mugy_io *ioman) {
  // Initialize the time independent factors.

  // Calculate the FLR operator building blocks: Gamma_0, <J_0>, hat{nabla}_\perp, hat{hat{nabla}}_\perp.
  struct mugy_flr *flr = mugy_flr_init(pop, grid);

  // Calculate the factors multiplying each moment in the Poisson equation.
  mugy_field_constop_init(pop, field, grid, flr);

  // Pre-compute the time independent linear operators in each equation.
  mugy_linear_constop_init(pop, grid, flr);

  // Precompute the FLR operators multipliying the potential inside Poisson brackets.
  mugy_nonlinear_constop_init(pop, grid, flr);

//    struct mugy_ad_file *fhr = mugy_io_create_population_perp_file(ioman, "pbFLRop", grid, pop, MUGY_REAL, MUGY_FOURIER_GRID, 3, 0);
//    mugy_io_write_mugy_array(NULL, "pbFLRop", fhr, pop->local->pbFLRop);
//    mugy_io_close_file(fhr);

//    struct mugy_ad_file *fhs = mugy_io_create_population_perp_file(ioman, "poissonDb", grid, pop, MUGY_REAL, MUGY_FOURIER_GRID, 1, 0);
//    mugy_io_write_mugy_array(NULL, "poissonDb", fhs, poissonDb);
//    mugy_io_close_file(fhs);
//    struct mugy_ad_file *fht = mugy_io_create_population_perp_file(ioman, "poissonSb", grid, pop, MUGY_REAL, MUGY_FOURIER_GRID, 1, 0);
//    mugy_io_write_mugy_array(NULL, "poissonSb", fht, poissonSb);
//    mugy_io_close_file(fht);

  // Free memory storing FLR building blocks (not needed anymore).
  mugy_flr_free(flr);
}

