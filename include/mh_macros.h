/* mugy: mh_macros.h
 *
 * A series of macros used throughout mugy.
 *
 */
#pragma once

#include "mh_userFLAGS.h"
#include <float.h> // For FLT_MIN, DBL_MIN.

// Number of dimensions in the code.
#define nDim 3

#ifdef USE_SINGLE_PRECISION
typedef float real;
#define MUGY_MPI_REAL MPI_FLOAT
#define MUGY_MPI_FOURIER MPI_C_COMPLEX
#define fmt_real "f"
static const real MUGY_REAL_MIN = FLT_MIN;
static const real MUGY_REAL_MAX = FLT_MAX;
#else
typedef double real;
#define MUGY_MPI_REAL MPI_DOUBLE
#define MUGY_MPI_FOURIER MPI_C_DOUBLE_COMPLEX
#define fmt_real "lf"
static const real MUGY_REAL_MIN = DBL_MIN;
static const real MUGY_REAL_MAX = DBL_MAX;
#endif

// Define our own int in case we wish to change to long.
typedef int mint;
#define MUGY_MPI_MINT MPI_INT
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
enum mugy_resource_mem {MUGY_HOST_MEM, MUGY_DEVICE_MEM, MUGY_HOSTDEVICE_MEM};
// Flags indicating whether to perform operation on host or device.
enum mugy_resource_calc {MUGY_HOST_CALC, MUGY_DEVICE_CALC};

// Types of data used in mugy.
enum mugy_data_types {MUGY_MINT, MUGY_REAL, MUGY_FOURIER};

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
enum mugy_memcpy_dir {
  MUGY_HOST2HOST     = cudaMemcpyHostToHost,
  MUGY_HOST2DEVICE   = cudaMemcpyHostToDevice,
  MUGY_DEVICE2HOST   = cudaMemcpyDeviceToHost,
  MUGY_DEVICE2DEVICE = cudaMemcpyDeviceToDevice,
};

#define DEFAULT_NUM_THREADS_DEV 256

#else

#define MUGY_CU_DH
#define MUGY_CU_D

// Directions of memcopy in host/device memory.
enum mugy_memcpy_dir {
  MUGY_HOST2HOST    ,
  MUGY_HOST2DEVICE  ,
  MUGY_DEVICE2HOST  ,
  MUGY_DEVICE2DEVICE,
};

#define DEFAULT_NUM_THREADS_DEV 1

#endif  // GPU related macros.
