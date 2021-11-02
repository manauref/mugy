/* mugy
   
   Functions used to allocate arrays.
*/
#include "alloc_mugy.h"

real* allocRealMoments(int ns, int nmom, struct gridType gridIn) {
  return (real *) calloc(ns*nmom*prod(gridIn.NxG), sizeof(real));
}

fourier *allocFourierMoments(int ns, int nmom, struct gridType gridIn) {
  return (fourier *) calloc(ns*nmom*prod(gridIn.NkxG), sizeof(fourier));
}

real* allocAliasedRealMoments(int ns, int nmom, struct gridType gridIn) {
  return (real *) calloc(ns*nmom*prod(gridIn.NxaG), sizeof(real));
}

fourier *allocAliasedFourierMoments(int ns, int nmom, struct gridType gridIn) {
  return (fourier *) calloc(ns*nmom*prod(gridIn.NkxaG), sizeof(fourier));
}
