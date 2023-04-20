/* mugy: reduce.c
 *
 * Functions to perform reductions.
 *
 */
#include "mh_reduce.h"

bool mugy_reduce_array_isfinite(struct mugy_array *arr) {
  // Check that all elements in a mugy_array are finite.
#ifdef USE_GPU
  enum mugy_resource_calc onResource = MUGY_DEVICE_CALC;
#else
  enum mugy_resource_calc onResource = MUGY_HOST_CALC;
#endif

  return mugy_array_isfinite(arr, onResource);
}
