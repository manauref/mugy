/* mugy: mh_field_dev.h
 *
 * Manipulate the field object on the device.
 *
 */
#pragma once

#include "mh_field.h"

// Solve the Poisson equation to obtain phik using the charge density
// in the time-stepping index 'tstepIdx'
void mugy_field_poisson_solve_dev(struct mugy_field *field, struct mugy_population *pop, struct mugy_grid *grid, mint tstepIdx);
