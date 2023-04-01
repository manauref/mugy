/* mugy: mh_alloc.h

   Functions used to allocate arrays.
*/

#ifndef ALLOC_MUGY
#define ALLOC_MUGY

#include "mh_data.h"
#include "mh_utilities.h"

// Wrappers to basic functions that allocate plain arrays.
//   numElements: number of elements in the array.
mint* alloc_mintArray(mint numElements);
char* alloc_charArray(mint numElements);
real* alloc_realArray(mint numElements);
fourier* alloc_fourierArray(mint numElements);

/* Allocate real-space moment vectors, on host and/or device.
     grid: grid on which to allocate the vector of moments.
     pop: population struct containing the number of species and moments.
     res: resource on which to allocate (host, device or both).
     mom: struct holding the vector of moments. */
void alloc_realMoments(const struct realGrid grid, const struct population pop, const resource res, struct realMoments *mom);

/* Allocate Fourier-space moment vectors, on host and/or device.
     grid: grid on which to allocate the vector of moments.
     pop: population struct containing the number of species and moments.
     res: resource on which to allocate (host, device or both).
     momk: struct holding the vector of moments. */
void alloc_fourierMoments(const struct fourierGrid grid, const struct population pop, const resource res, struct fourierMoments *momk);

/* Functions to free memory associated used for moment vector. 
     mom: vector of moments to free.
     res: resource where memory needs to be freed (host, device or both). */
void free_realMoments(struct realMoments *mom, const resource res);
void free_fourierMoments(struct fourierMoments *momk, const resource res);

#endif
