/* A mugy header file:

   Data types (e.g. structs) used in mugy.
*/
#ifndef DATA_MUGY
#define DATA_MUGY

#include <complex.h>  /* For complex data types. */
#include <stdbool.h>  // e.g. for bool, true, false.

#if USE_SINGLE_PRECISION > 0
typedef float real;
typedef float complex fourier;
#else
typedef double real;
typedef double complex fourier;
#endif

// Container for IO instructions
struct ioSetup {
  char *inputFile;           // Name of input file.
  char *outputDir;           // Address of output directory.
  char *restartDir;          // Address of restart directory.
  bool isRestart;            // Is this simulation a restart of a previous one?
  bool outToOldDir;          // If restart, output to directory of previous run?
};

// Flags for differentiating between operations done on device
// or on host, or both.
typedef enum {hostOnly, deviceOnly, hostAndDevice} resource;

struct gridType {
  int NxaG[3];    // Number of cells in aliased configuration space.
  int NkxaG[3];   // Number of distinct (absolute) aliased wavenumbers.
  int NekxaG[3];  // Number of elements in an aliased k-space array.
  int NkxG[3];    // Number of distinct absolute amplitude wavenumbers (de-aliased).
  int NekxG[3];   // Number of elements in a de-aliased k-space array.
  int NxG[3];     // Number of cells in de-aliased configuration space.
  real kxMin[3];  // Minimum finite absolute amplitude wavenumbers.
};

// Structures storing host and device pointers to vector field.
struct realMoments {
  real *ho;
  real *dev;
};
struct fourierMoments {
  fourier *ho;
  fourier *dev;
};

#endif
