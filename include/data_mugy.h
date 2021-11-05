/* A mugy header file:

   Data types (e.g. structs) used in mugy.
*/
#ifndef DATA_MUGY
#define DATA_MUGY

#include <complex.h>  /* For complex data types. */
#include <stdbool.h>  // e.g. for bool, true, false.
#include "mugyFLAGS.h"

// Number of dimensions in the code.
#define nDim 3

#if USE_SINGLE_PRECISION > 0
typedef float real;
typedef float complex fourier;
#define SCNfREAL "f"
#else
typedef double real;
typedef double complex fourier;
#define SCNfREAL "lf"
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

struct realGrid {
  int Nx[3];     // Number of cells.
  real xMin[3];  // Minimum coordinates.
  real xMax[3];  // Maximum coordinates.
  real dx[3];    // Cell length.
  real *x;       // Coordinates in each direction.
};

struct fourierGrid {
  int Nkx[3];           // Number of distinct absolute amplitude wavenumbers (counting k=0).
  int Nekx[3];          // Number of elements in a k-space array (counting k=0 and negative k's).
  real kxMin[3];        // Minimum finite absolute amplitude wavenumbers.
  real *kx;             // Coordinates along each direction.
  struct realGrid dual; // Real grid dual to this Fourier grid.
  real kxMaxDyn[3];     // Maximum k evolved. Above this we multiply time rates of change by zero.
};

struct grid {
  struct fourierGrid fG;  // this grid's (dealised) Fourier grid.
  struct fourierGrid fGa; // this grid's aliased Fourier grid.
  int mpiProcs[3];        // Number of MPI processes along each direction.
};

struct timeSetup {
  real dt;               // Time step.
  real endTime;          // Absolute end time (from t=0 of first simulation).
  int nFrames;           // Absolute frames to output (from t=0 of first simulation).
  int ark_kySplit;       // multirate splitting index: highest "slow" wavenumber (higher are "fast").
  int ark_fastTableExp;  // Butcher table index for fast explicit method.
  int ark_fastTableImp;  // Butcher table index for fast implicit method.
  int ark_slowTable;     // Butcher table index for slow method.
  real ark_dtFast;       // fixed 'fast' time step to use.
  real ark_rtol;         // relative solution tolerance for temporal adaptivity.
  real ark_atol;         // absolute solution tolerance for temporal adaptivity
  int ark_ewtScaling;    // which error weight scaling to use.
};

struct speciesParameters {
  int numSpecies;    // Number of species.
  int numMoments;    // Number of moments.
  int mpiProcs;      // MPI decomposition of the species.
  real *qCharge;     // Charge.
  real *muMass;      // sqrt of the mass.
  real *tau;         // Temperature. 
  real *omSt;        // omega_star.
  real *omd;         // omega_d.
  real *delta;       // (Tpar + Tperp)/T.
  real *deltaPerp;   // Tperp/T.
  real *eta;         // L_n/L_{Tperp}.
  real *alpha;       // Damping coefficient.
  real *nu;          // Diffusion coefficient.
  real *delta0;      // i-delta non-adiabaticity parameter.
  real *hDiffOrder;  // Hyperdiffusion order.
  real *hDiff;       // Hyperdiffusion coefficient.
  real *kDiffMin;    // Minimum k at which to apply HD.
  // The following are used by initial conditions.
  int *icOp;         // IC option.
  real *initAux;  // Auxiliary parameters for ICs.
  real *initA;       // Initial amplitude.
  real *noiseA;      // Initial noise amplitude.
};

struct fieldParameters {
  real lambdaD;  // Debye shielding parameter (normalized Debye length).
  int pade;      // Option to use Pade approximations.
  // The following are used by initial conditions.
  int icOp;      // IC option.
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

// Define the various fields and moments needed.
extern struct realMoments mom, moma;
extern struct fourierMoments momk, momka;

#endif
