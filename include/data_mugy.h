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
#define mpi_real MPI_FLOAT
#define mpi_fourier MPI_C_COMPLEX
#define fmt_real "f"
#else
typedef double real;
typedef double complex fourier;
#define fmt_real "lf"
#define mpi_real MPI_DOUBLE
#define mpi_fourier MPI_C_DOUBLE_COMPLEX
#endif

// Define our own int in case we wish to change to long.
typedef int mint;
#define mpi_mint MPI_INT
#define fmt_mint "d"

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
  mint Nx[nDim];        // Number of cells.
  mint NxTot;           // Total number of cells.
  mint NxyTot;          // Total number of cells in an x-y plane.
  real xMin[nDim];      // Minimum coordinates.
  real xMax[nDim];      // Maximum coordinates.
  real dx[nDim];        // Cell length.
  real *x;              // Coordinates in each direction.
  mint globalOff[nDim]; // Offset of first element in this process within the global domain.
};

struct fourierGrid {
  mint Nkx[nDim];        // Number of distinct absolute amplitude wavenumbers (counting k=0).
  mint Nekx[nDim];       // Number of elements in a k-space array (counting k=0 and negative k's).
  mint NekxTot;          // Total number of elements.
  mint NekxyTot;         // Total number of cells in an kx-ky plane.
  real kxMin[nDim];      // Minimum finite absolute amplitude wavenumbers.
  real *kx;              // Coordinates along each direction.
  struct realGrid dual;  // Real grid dual to this Fourier grid.
  real kxMaxDyn[nDim];   // Maximum k evolved. Above this we multiply time rates of change by zero.
  mint globalOff[nDim];  // Offset of first element in this process within the global domain.
};

struct grid {
  struct fourierGrid fG;  // This grid's (dealised) Fourier grid.
  struct fourierGrid fGa; // This grid's aliased Fourier grid.
  mint mpiProcs[nDim];    // Number of MPI processes along each direction.
};

struct timeSetup {
  real dt;               // Time step.
  real endTime;          // Absolute end time (from t=0 of first simulation).
  mint nFrames;          // Absolute frames to output (from t=0 of first simulation).
  mint ark_kySplit;      // multirate splitting index: highest "slow" wavenumber (higher are "fast").
  mint ark_fastTableExp; // Butcher table index for fast explicit method.
  mint ark_fastTableImp; // Butcher table index for fast implicit method.
  mint ark_slowTable;    // Butcher table index for slow method.
  real ark_dtFast;       // fixed 'fast' time step to use.
  real ark_rtol;         // relative solution tolerance for temporal adaptivity.
  real ark_atol;         // absolute solution tolerance for temporal adaptivity
  mint ark_ewtScaling;   // which error weight scaling to use.
};

struct species {
  mint numMoments;   // Number of moments.
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
  mint icOp;         // IC option.
  real *initAux;     // Auxiliary parameters for ICs.
  real initA;        // Initial amplitude.
  real noiseA;       // Initial noise amplitude.
};

struct population {
  mint numSpecies;       // Number of species.
  mint mpiProcs;         // MPI decomposition of the species.
  mint globalSpecOff;    // Offset of first species in this process within the global population.
  mint globalMomOff;     // Offset of first moment in this process within global number of moments.
  struct species *spec;  // Pointer to array of species.
  mint numMomentsTot;    // Total number of moments across all species.
};

struct fieldParameters {
  real lambdaD;  // Debye shielding parameter (normalized Debye length).
  mint pade;      // Option to use Pade approximations.
  // The following are used by initial conditions.
  mint icOp;      // IC option.
};

struct timeState {
  real simTime;
  mint time;
  mint framesOut;
  mint hdAdjusts;
  mint dtAdjusts;
  real dt;
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
fourier* getMoment_fourier(struct fourierGrid grid, struct population pop, const mint sIdx, const mint momIdx, fourier *momkIn);

// Linear index given the nDim-dimensional subscript in a Fourier grid.
mint sub2lin_fourier(const mint *kxI, const struct fourierGrid grid);

// nDim-dimensional subscript given the linear index in a Fourier grid.
void lin2sub_fourier(mint *kxI, mint lin, const struct fourierGrid grid);

// (kx,ky,kz) coordinates given the multidimensional kxI index.
void get_kx(real *kx, mint *kxI, const struct fourierGrid grid);

#endif
