/* A mugy header file:

   Functions used to allocate arrays.
*/

#ifndef ALLOC_MUGY
#define ALLOC_MUGY

#include "data_mugy.h"
#include "utilities.h"
#include <stdlib.h>  // e.g. for calloc.

// Wrappers to basic functions that allocate plain arrays.
//   numElements: number of elements in the array.
real* alloc_realArray(int numElements);
fourier *alloc_fourierArray(int numElements);

/* Allocate real-space moment vectors, on host and/or device.
     ns: number of species.
     nmom: number of moments per species.
     grid: grid on which to allocate the vector of moments.
     res: resource on which to allocate (host, device or both).
     mom: struct holding the vector of moments. */
void alloc_realMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct realMoments *mom);

/* Allocate Fourier-space moment vectors, on host and/or device.
     ns: number of species.
     nmom: number of moments per species.
     grid: grid on which to allocate the vector of moments.
     res: resource on which to allocate (host, device or both).
     momk: struct holding the vector of moments. */
void alloc_fourierMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct fourierMoments *momk);

/* Allocate aliased real-space moment vectors, on host and/or device.
     ns: number of species.
     nmom: number of moments per species.
     grid: grid on which to allocate the vector of moments.
     res: resource on which to allocate (host, device or both).
     mom: struct holding the vector of moments. */
void alloc_aliasedRealMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct realMoments *mom);

/* Allocate aliased Fourier-space moment vectors, on host and/or device.
     ns: number of species.
     nmom: number of moments per species.
     grid: grid on which to allocate the vector of moments.
     res: resource on which to allocate (host, device or both).
     momk: struct holding the vector of moments. */
void alloc_aliasedFourierMoments(const int ns, const int nmom, const struct gridType grid, const resource res, struct fourierMoments *momk);


/* Functions to free memory associated used for moment vector. 
     mom: vector of moments to free.
     res: resource where memory needs to be freed (host, device or both). */
void free_realMoments(struct realMoments *mom, const resource res);
void free_fourierMoments(struct fourierMoments *momk, const resource res);

#endif
