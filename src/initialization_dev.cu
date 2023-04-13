/* mugy: initialization_dev.c
 *
 * Functions used to initialize the device (GPU).
 *
 */

extern "C" {
#include "mh_initialization_dev.h"
#include "mh_utilities_dev.h"
}

void device_init_dev(struct mugy_comms *comms) {

  // Set the device of this MPI rank.
  int devCount = 0;
  checkCudaErrors(cudaGetDeviceCount(&devCount));
  mint mpiRank = comms->world.rank;
  mint devID = mpiRank % devCount;
  checkCudaErrors(cudaSetDevice(devID));
  printf("  My rank: %d | # of GPUs: %d | my GPU: %d\n", mpiRank, devCount, devID);
  printf("\n");

  // Get device properties.
  cudaDeviceProp deviceProp;
  checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

  printf("  Device properties:\n");
  printf("    global memory = %zu\n",deviceProp.totalGlobalMem);
  printf("    shared memory/block = %zu\n",deviceProp.sharedMemPerBlock);
  printf("    warp size = %d\n",deviceProp.warpSize);
  printf("    max threads/block = %d\n",deviceProp.maxThreadsPerBlock);
  printf("    major compute capability = %d\n",deviceProp.major);
  printf("    minor compute capability = %d\n",deviceProp.minor);
  printf("\n");

}

