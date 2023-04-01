/* mugy: mh_macros.h

   A series of macros used throughout mugy.
*/

#ifndef MUGY_MACROS 
#define MUGY_MACROS 


// GPU related macros. USE_GPU is passed as a compile-time preprocessor variable.
#ifdef USE_GPU

#include <cuda_runtime.h>

#define CU_DH __device__ __host__
#define CU_D __device__

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host = cudaMemcpyHostToHost,
  host2dev  = cudaMemcpyHostToDevice,
  dev2host  = cudaMemcpyDeviceToHost,
  dev2dev   = cudaMemcpyDeviceToDevice,
};

#define DEFAULT_NUM_DEV_THREADS 256

#else

#define CU_DH
#define CU_D

// Directions of memcopy in host/device memory.
enum memcpy_dir_dev {
  host2host,
  host2dev ,
  dev2host ,
  dev2dev  ,
};

#define DEFAULT_NUM_DEV_THREADS 1

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


#endif
