/* mugy: mh_macros.h

   A series of macros used throughout mugy.
*/

#ifndef MUGY_MACROS 
#define MUGY_MACROS 

// Number of dimensions in the code.
#define nDim 3

#if USE_SINGLE_PRECISION > 0
typedef float real;
#else
typedef double real;
#endif

// Define our own int in case we wish to change to long.
typedef int mint;

// Moment indices.
#define denIdx 0
#define tempIdx 1

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


// GPU related macros. USE_GPU is passed as a compile-time preprocessor variable.
#if USE_GPU

#include <cuda_runtime.h>

#define CU_DH __device__ __host__
#define CU_D __device__

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host     = cudaMemcpyHostToHost,
  host2device   = cudaMemcpyHostToDevice,
  device2host   = cudaMemcpyDeviceToHost,
  device2device = cudaMemcpyDeviceToDevice,
};

#define DEFAULT_NUM_THREADS_DEV 256

#else

#define CU_DH
#define CU_D

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host    ,
  host2device  ,
  device2host  ,
  device2device,
};

#define DEFAULT_NUM_THREADS_DEV 1

#endif  // GPU related macros.


#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

// Number of fields needed for time stepping.
#if (TIME_STEPPER == 4)
#define TIME_STEPPER_NUM_FIELDS 4
#endif

#endif
