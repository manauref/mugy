/* mugy: mh_macros.h

   A series of macros used throughout mugy.
*/

#ifndef MUGY_MACROS 
#define MUGY_MACROS 

#include "mh_userFLAGS.h"

// Number of dimensions in the code.
#define nDim 3

#if USE_SINGLE_PRECISION
typedef float real;
#define mpi_real MPI_FLOAT
#define mpi_fourier MPI_C_COMPLEX
#define fmt_real "f"
#else
typedef double real;
#define mpi_real MPI_DOUBLE
#define mpi_fourier MPI_C_DOUBLE_COMPLEX
#define fmt_real "lf"
#endif

// Define our own int in case we wish to change to long.
typedef int mint;
#define mpi_mint MPI_INT
#define fmt_mint "d"

// Moment indices.
#define denIdx 0
#define tempIdx 1

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// ID of rank that does simpe I/O:
#define ioRank 0

// Flag indicating whether to use host or device memory, or both.
enum resource_mem {hostMem, deviceMem, hostAndDeviceMem};
// Flags indicating whether to perform operation on host or device.
enum resource_comp {defaultComp, hostComp, deviceComp};

// Types of data used in mugy.
enum mugy_datatype {mint_enum, real_enum, fourier_enum};

#define mugy_max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define mugy_min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

// Number of fields needed for time stepping.
#if (TIME_STEPPER == 4)
#define TIME_STEPPER_NUM_FIELDS 4
#endif


// GPU related macros. USE_GPU is passed as a compile-time preprocessor variable.
#if USE_GPU

#include <cuda_runtime.h>

#define MUGY_CU_DH __device__ __host__
#define MUGY_CU_D __device__

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host     = cudaMemcpyHostToHost,
  host2device   = cudaMemcpyHostToDevice,
  device2host   = cudaMemcpyDeviceToHost,
  device2device = cudaMemcpyDeviceToDevice,
};

#define DEFAULT_NUM_THREADS_DEV 256

#else

#define MUGY_CU_DH
#define MUGY_CU_D

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host    ,
  host2device  ,
  device2host  ,
  device2device,
};

#define DEFAULT_NUM_THREADS_DEV 1

#endif  // GPU related macros.


#endif
