/* A mugy header file:

   Functions used to allocate arrays.
*/

#ifndef ALLOC_MUGY
#define ALLOC_MUGY

#include "data_mugy.h"
#include "utilities.h"
#include <stdlib.h>  // e.g. for calloc.

/* Allocate a vector of real-space moments.
     ns: number of species.
     nmom: number of moments per species.
     grodIn: grid on which to allocate the vector of moments. */
real* allocRealMoments(int ns, int nmom, struct gridType gridIn);

/* Allocate a vector of Fourier-space moments.
     ns: number of species.
     nmom: number of moments per species.
     grodIn: grid on which to allocate the vector of moments. */
fourier *allocFourierMoments(int ns, int nmom, struct gridType gridIn);

/* Allocate a vector of aliased real-space moments.
     ns: number of species.
     nmom: number of moments per species.
     grodIn: grid on which to allocate the vector of moments. */
real* allocAliasedRealMoments(int ns, int nmom, struct gridType gridIn);

/* Allocate a vector of aliased Fourier-space moments.
     ns: number of species.
     nmom: number of moments per species.
     grodIn: grid on which to allocate the vector of moments. */
fourier *allocAliasedFourierMoments(int ns, int nmom, struct gridType gridIn);

#endif
