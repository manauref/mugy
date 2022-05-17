/* mugy: macros.h

   A series of macros used throughout mugy.
*/


// GPU related macros. USE_GPU is passed as a compile-time preprocessor variable.
#ifdef USE_GPU

#include <cuda_runtime.h>

#define CU_DH __device__ __host__
#define CU_D __device__

// for directional copies
enum dev_memcpy_dir {
  DEV_MEMCPY_H2D = cudaMemcpyHostToDevice,
  DEV_MEMCPY_D2H = cudaMemcpyDeviceToHost,
  DEV_MEMCPY_D2D = cudaMemcpyDeviceToDevice
};

#define DEFAULT_NUM_DEV_THREADS 256

#else

#define CU_DH
#define CU_D

// for directional copies
enum dev_memcpy_dir {
  DEV_MEMCPY_H2D,
  DEV_MEMCPY_D2H,
  DEV_MEMCPY_D2D,
};

#define DEFAULT_NUM_DEV_THREADS 1

#endif // CUDA specific defines etc


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
