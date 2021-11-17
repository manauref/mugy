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
#define mugyMPI_REAL MPI_FLOAT
#define mugyMPI_FOURIER MPI_C_COMPLEX
#define SCNfREAL "f"
#else
typedef double real;
typedef double complex fourier;
#define SCNfREAL "lf"
#define mugyMPI_REAL MPI_DOUBLE
#define mugyMPI_FOURIER MPI_C_DOUBLE_COMPLEX
#endif

#define mugyMPI_INT MPI_INT

// ID of rank that does simpe I/O:
#define ioRank 0

// Moment indices.
#define denIdx 0 
#define tempIdx 1

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
  int Nx[nDim];         // Number of cells.
  int NxTot;            // Total number of cells.
  int NxyTot;           // Total number of cells in an x-y plane.
  real xMin[nDim];      // Minimum coordinates.
  real xMax[nDim];      // Maximum coordinates.
  real dx[nDim];        // Cell length.
  real *x;              // Coordinates in each direction.
  int globalOff[nDim];  // Offset of first element in this process within the global domain.
};

struct fourierGrid {
  int Nkx[nDim];         // Number of distinct absolute amplitude wavenumbers (counting k=0).
  int Nekx[nDim];        // Number of elements in a k-space array (counting k=0 and negative k's).
  int NekxTot;           // Total number of elements.
  int NekxyTot;          // Total number of cells in an kx-ky plane.
  real kxMin[nDim];      // Minimum finite absolute amplitude wavenumbers.
  real *kx;              // Coordinates along each direction.
  struct realGrid dual;  // Real grid dual to this Fourier grid.
  real kxMaxDyn[nDim];   // Maximum k evolved. Above this we multiply time rates of change by zero.
  int globalOff[nDim];   // Offset of first element in this process within the global domain.
};

struct grid {
  struct fourierGrid fG;  // This grid's (dealised) Fourier grid.
  struct fourierGrid fGa; // This grid's aliased Fourier grid.
  int mpiProcs[nDim];     // Number of MPI processes along each direction.
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

struct species {
  int numMoments;    // Number of moments.
  real qCharge;      // Charge.
  real muMass;       // sqrt of the mass.
  real tau;          // Temperature. 
  real omSt;         // omega_star.
  real omd;          // omega_d.
  real delta;        // (Tpar + Tperp)/T.
  real deltaPerp;    // Tperp/T.
  real eta;          // L_n/L_{Tperp}.
  real *alpha;       // Damping coefficient.
  real *nu;          // Diffusion coefficient.
  real delta0;       // i-delta non-adiabaticity parameter.
  real hDiffOrder;   // Hyperdiffusion order.
  real *hDiff;       // Hyperdiffusion coefficient.
  real *kDiffMin;    // Minimum k at which to apply HD.
  // The following are used by initial conditions.
  int icOp;          // IC option.
  real *initAux;     // Auxiliary parameters for ICs.
  real initA;        // Initial amplitude.
  real noiseA;      // Initial noise amplitude.
};

struct population {
  int numSpecies;        // Number of species.
  int mpiProcs;          // MPI decomposition of the species.
  int globalOff;         // Offset of first element in this process within the global domain.
  struct species *spec;  // Pointer to array of species.
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

// Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
fourier* getMoment_fourier(struct fourierGrid grid, struct population pop, const int sIdx, const int momIdx, fourier *momkIn);

// Linear index given the nDim-dimensional subscript in a Fourier grid.
int sub2lin_fourier(const int *kxI, const struct fourierGrid grid);

// nDim-dimensional subscript given the linear index in a Fourier grid.
void lin2sub_fourier(int *kxI, int lin, const struct fourierGrid grid);

// (kx,ky,kz) coordinates given the multidimensional kxI index.
void get_kx(real *kx, int *kxI, const struct fourierGrid grid);

#endif
