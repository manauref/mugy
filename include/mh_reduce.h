/* mugy: mh_reduce.h
 *
 * Functions to perform reductions.
 *
 */
#pragma once

#include "mh_array.h"

// Check that all elements in a mugy_array are finite.
bool mugy_reduce_array_isfinite(struct mugy_array *arr);
