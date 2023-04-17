/* mugy: time.c
 *
 * Functions to manipulate the time struct and do time integration.
 *
 */
#include "mh_time.h"
#include <stdlib.h>  // for malloc.

struct mugy_time *mugy_time_alloc() {
  // Allocate the time object.
  struct mugy_time *time = (struct mugy_time *) malloc(sizeof(struct mugy_time));
  return time;
}

void mugy_time_free(struct mugy_time *time) {
  free(time);
}
