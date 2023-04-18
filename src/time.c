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

void mugy_time_wcstamp(struct timespec *wcstamp) {
  // Record the wall-clock time in a wc stamp.
  clock_gettime(CLOCK_MONOTONIC, wcstamp);
};

double mugy_time_sec(struct timespec wcstamp) {
  // Return the time of a wall-clock time in seconds.
  return wcstamp.tv_sec + 1e-9*wcstamp.tv_nsec;
}

struct timespec mugy_time_diff(struct timespec start, struct timespec end) {
  // Compute the time difference between two wall-clock stamps.
  struct timespec tm;
  if ((end.tv_nsec - start.tv_nsec)<0) {
    tm.tv_sec  = end.tv_sec-start.tv_sec-1;
    tm.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    tm.tv_sec  = end.tv_sec-start.tv_sec;
    tm.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return tm;
}

double mugy_time_elapsed_sec(struct timespec prev) {
  // Return the time elapsed since the a previous wall-clock stamp (prev).
  struct timespec curr;
  mugy_time_wcstamp(&curr);
  return mugy_time_sec(mugy_time_diff(prev, curr));
}

void mugy_time_free(struct mugy_time *time) {
  // Free memory associated with time object.
  free(time);
}
