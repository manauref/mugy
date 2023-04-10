/* mugy: mh_field.h
 *
 * Field object holding the electromagnetic fields (e.g. phi)
 * and functions to use/change it.
 *
 */
#ifndef MUGY_FIELD
#define MUGY_FIELD

#include "mh_data.h"

struct mugy_fieldParameters {
  real lambdaD;  // Debye shielding parameter (normalized Debye length).
  mint pade;      // Option to use Pade approximations.
  // The following are used by initial conditions.
  mint icOp;      // IC option.
};

#endif
