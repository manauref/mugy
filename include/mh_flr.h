/* mugy: mh_flr.h
 *
 * Create various FLR (finite Larmor radius) operators.
 *
 */
#pragma once

#include "mh_population.h"
#include "mh_grid.h"
#include "mh_field.h"
#include "mh_io_adios.h"

void mugy_flr_init(struct mugy_population *pop, struct mugy_grid *grids, struct mugy_field *field, struct mugy_ioManager *ioman);
