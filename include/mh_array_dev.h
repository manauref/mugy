/* mugy: mh_array_dev.h
 *  
 * Device methods for a mugy_array.
 *
 */
#pragma once

#include "mh_array.h"

// Set all elements in the array to zero.
void mugy_array_zero_dev(struct mugy_array *arr);

// Scale an array by a factor 'fac'.
void mugy_array_scale_dev(struct mugy_array *arr, real fac);

// Check that none of the elements are inf or NaN.
bool mugy_array_isfinite_dev(struct mugy_array *arr);
