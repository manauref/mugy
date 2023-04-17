/* mugy: constop.c
 *
 * Precompute time independent operators.
 *
 */
#pragma once

#include "mh_population.h"
#include "mh_grid.h"
#include "mh_field.h"
#include "mh_io_adios.h"

// Initialize the time independent factors.
void mugy_constop_init(struct mugy_population *pop, struct mugy_grid *grid, struct mugy_field *field, struct mugy_io *ioman);
