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

// Flag indicating whether to use host or device memory, or both.
enum resource_mem {hostMem, deviceMem, hostAndDeviceMem};
// Flags indicating whether to perform operation on host or device.
enum resource_comp {defaultComp, hostComp, deviceComp};

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
