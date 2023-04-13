/* mugy: mh_field.h
 *
 * Field object holding the electromagnetic fields (e.g. phi)
 * and functions to use/change it.
 *
 */
#pragma once

#include "mh_macros.h"
#include "mh_array.h"
#include "mh_grid.h"
#include "mh_population.h"

struct mugy_field_pars {
  real lambdaD;  // Debye shielding parameter (normalized Debye length).
  mint pade;     // Option to use Pade approximations.
  // The following are used by initial conditions.
  mint icOp;     // IC option.
};

struct mugy_field {
  struct mugy_field_pars pars;  // Field parameters.
  struct mugy_array *phik;      // Potential in Fourier space.
  struct mugy_array *gyrophik;  // Gyroaveraged potentials in Fourier space.
};

// Allocate the field object.
struct mugy_field *mugy_field_alloc();

// Initialize the rest of the field arrays.
void mugy_field_init(struct mugy_field *field, struct mugy_grid *grid, struct mugy_population *pop);

// Free the field object.
void mugy_field_free(struct mugy_field *field);
