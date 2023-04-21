/* mugy: mh_array_dev.h
 *  
 * Device methods for a mugy_array.
 *
 */
#pragma once

#include "mh_array.h"

// Set all elements in the array to zero.
void mugy_array_zero_dev(struct mugy_array *arr);

// Increment array out by array x.
void mugy_array_increment_dev(struct mugy_array *out, struct mugy_array *x);

// Assign array out with a*x.
void mugy_array_ax_assign_dev(struct mugy_array *out, real a, struct mugy_array *x);

// Assign array out with a*x+y.
void mugy_array_axpy_assign_dev(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y);

// Increment array out with a*x+y.
void mugy_array_axpy_increment_dev(struct mugy_array *out, real a, struct mugy_array *x, struct mugy_array *y);

// Scale an array by a factor 'fac'.
void mugy_array_scale_dev(struct mugy_array *arr, real fac);

// Check that none of the elements are inf or NaN.
bool mugy_array_isfinite_dev(struct mugy_array *arr);
