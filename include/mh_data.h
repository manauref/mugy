/* mugy: mh_data.h

   Data types (e.g. structs) used in mugy.
*/
#ifndef MUGY_DATA
#define MUGY_DATA

#include <complex.h>  /* For complex data types. */
#include <stdbool.h>  // e.g. for bool, true, false.
#include "mh_userFLAGS.h"
#include "mh_macros.h"

// mugy specific types (some other types in mh_macros.h).
#if USE_SINGLE_PRECISION
typedef float complex fourier;
#define mpi_real MPI_FLOAT
#define mpi_fourier MPI_C_COMPLEX
#define fmt_real "f"
#else
typedef double complex fourier;
#define fmt_real "lf"
#define mpi_real MPI_DOUBLE
#define mpi_fourier MPI_C_DOUBLE_COMPLEX
#endif

#define mpi_mint MPI_INT
#define fmt_mint "d"

// Container for IO instructions
struct mugy_ioSetup {
  char *inputFile;           // Name of input file.
  char *outputDir;           // Address of output directory.
  char *restartDir;          // Address of restart directory.
  bool isRestart;            // Is this simulation a restart of a previous one?
  bool outToOldDir;          // If restart, output to directory of previous run?
};

// Flag indicating whether to use host or device memory, or both.
enum resource_mem {hostMem, deviceMem, hostAndDeviceMem};
// Flags indicating whether to perform operation on host or device.
enum resource_comp {defaultComp, hostComp, deviceComp};

struct mugy_realGrid {
  mint Nx[nDim];        // Number of cells.
  mint NxTot;           // Total number of cells.
  mint NxyTot;          // Total number of cells in an x-y plane.
  real xMin[nDim];      // Minimum coordinates.
  real xMax[nDim];      // Maximum coordinates.
  real dx[nDim];        // Cell length.
  real *x;              // Coordinates in each direction.
  mint globalOff[nDim]; // Offset of first element in this process within the global domain.
};

struct mugy_fourierGrid {
  mint Nkx[nDim];        // Number of distinct absolute amplitude wavenumbers (counting k=0).
  mint Nekx[nDim];       // Number of elements in a k-space array (counting k=0 and negative k's).
  mint NekxTot;          // Total number of elements.
  mint NekxyTot;         // Total number of cells in an kx-ky plane.
  real kxMin[nDim];      // Minimum finite absolute amplitude wavenumbers.
  real *kx;              // Coordinates along each direction.
  struct mugy_realGrid dual;  // Real grid dual to this Fourier grid.
  real kxMaxDyn[nDim];   // Maximum k evolved. Above this we multiply time rates of change by zero.
  mint globalOff[nDim];  // Offset of first element in this process within the global domain.
};

struct mugy_grid {
  struct mugy_fourierGrid fG;  // This grid's (dealised) Fourier grid.
  struct mugy_fourierGrid fGa; // This grid's aliased Fourier grid.
  mint mpiProcs[nDim];    // Number of MPI processes along each direction.
};

struct mugy_timeSetup {
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

struct mugy_species {
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
  real *hDiffOrder;  // Hyperdiffusion order.
  real *hDiff;       // Hyperdiffusion coefficient.
  real *kDiffMin;    // Minimum k at which to apply HD.
  // The following are used by initial conditions. It'd be good to put them elsewhere.
  mint icOp;         // IC option.
  real *initAux;     // Auxiliary parameters for ICs.
  real initA;        // Initial amplitude.
  real noiseA;       // Initial noise amplitude.
};

// Structures storing an array on host, device, or both.
struct mugy_realArray {
  real *ho;    // Pointer to host memory.
  real *dev;   // Pointer to device memory.
  mint nelem;  // Number of elements allocated.
};
struct mugy_fourierArray {
  fourier *ho;   // Pointer to host memory.
  fourier *dev;  // Pointer to device memory.
  mint nelem;    // Number of elements allocated.
};

struct mugy_population {
  mint numSpecies;       // Number of species.
  mint mpiProcs;         // MPI decomposition of the species.
  mint globalSpecOff;    // Offset of first species in this process within the global population.
  mint globalMomOff;     // Offset of first moment in this process within global number of moments.
  struct mugy_species *spec;  // Pointer to array of species.
  mint numMomentsTot;    // Total number of moments across all species.
  struct mugy_fourierArray *momk; // Moments in Fourier space. Possibly multiple copies (e.g. for time stepper).
};

struct mugy_fieldParameters {
  real lambdaD;  // Debye shielding parameter (normalized Debye length).
  mint pade;      // Option to use Pade approximations.
  // The following are used by initial conditions.
  mint icOp;      // IC option.
};

struct mugy_timeState {
  real simTime;
  mint time;
  mint framesOut;
  mint hdAdjusts;
  mint dtAdjusts;
  real dt;
};

// Return a pointer to the momIdx-th moment of the sIdx-th species in mom/momk.
real* getMoment_real(struct mugy_realGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, real *momIn);
fourier* getMoment_fourier(struct mugy_fourierGrid grid, struct mugy_population pop, mint sIdx, mint momIdx, fourier *momkIn);

// Linear index given the nDim-dimensional subscript in a real/Fourier grid.
mint sub2lin_real(mint *xI, const struct mugy_realGrid grid);
mint sub2lin_fourier(mint *kxI, const struct mugy_fourierGrid grid);

// nDim-dimensional subscript given the linear index in a real/Fourier grid.
void lin2sub_real(mint *xI, mint lin, const struct mugy_realGrid grid);
void lin2sub_fourier(mint *kxI, mint lin, const struct mugy_fourierGrid grid);

// (x,y,z)/(kx,ky,kz) coordinates given the multidimensional xI/kxI index.
void get_x(real *x, mint *xI, const struct mugy_realGrid grid);
void get_kx(real *kx, mint *kxI, const struct mugy_fourierGrid grid);

// Copy mint/real/fourier data (between host and device, or within a host or device).
void memcpy_mint(mint *dest, mint *src, mint numElements, enum memcpy_dir_dev dir);
void memcpy_real(real *dest, real *src, mint numElements, enum memcpy_dir_dev dir);
void memcpy_fourier(void *dest, void *src, mint numElements, enum memcpy_dir_dev dir);

// Copy real/fourier array between host and device.
void hodevXfer_realArray(struct mugy_realArray *arr, enum memcpy_dir_dev dir);
void hodevXfer_fourierArray(struct mugy_fourierArray *arr, enum memcpy_dir_dev dir);

// Scale an array by a factor 'fac'.
void scale_realArray(struct mugy_realArray *arr, real fac, enum resource_comp res);
void scale_fourierArray(struct mugy_fourierArray *arrk, real fac, enum resource_comp res);

#endif
