/* mugy: data_mugy
   
   Parameters and fields used throughout mugy.
*/

#include "mh_data.h"
#include "mh_data_dev.h"
#include <string.h>

real* getMoment_real(struct realGrid grid, struct population pop, mint sIdx, mint momIdx, real *momIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in mom.
  real* ptrOut = momIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NxTot;
}

fourier* getMoment_fourier(struct fourierGrid grid, struct population pop, mint sIdx, mint momIdx, fourier *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NekxTot;
}

mint sub2lin_real(mint *xI, const struct realGrid grid) {
  // Given the nDim-dimensional index (subscript) xI return the linear index
  // in a real grid. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid.Nx[1],1,grid.NxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += xI[d]*strides[d];
  return lin;
}

mint sub2lin_fourier(mint *kxI, const struct fourierGrid grid) {
  // Given the nDim-dimensional index (subscript) kxI return the linear index
  // in a Fourier grid. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += kxI[d]*strides[d];
  return lin;
}

void lin2sub_real(mint *xI, mint lin, const struct realGrid grid) {
  // Given the linear index 'lin' in a real grid, return the nDim-dimensional
  // index (subscript) xI. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid.Nx[1],1,grid.NxyTot};
  for (mint d=0; d<nDim; d++) {
    xI[d] = lin/strides[d];
    lin -= xI[d]*strides[d];
  }
}

void lin2sub_fourier(mint *kxI, mint lin, const struct fourierGrid grid) {
  // Given the linear index 'lin' in a Fourier grid, return the nDim-dimensional
  // index (subscript) kxI. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  for (mint d=0; d<nDim; d++) {
    kxI[d] = lin/strides[d];
    lin -= kxI[d]*strides[d];
  }
}

void get_x(real *x, mint *xI, const struct realGrid grid) {
  // Obtain the x=(x,y,z) coordinates given the multidimensional
  // xI index. Assume the flat grid.x array is organized as {x,y,z}.
  real* x_p = grid.x;
  for (mint d=0; d<nDim; d++) {
    x[d] = x_p[xI[d]]; 
    x_p += grid.Nx[d];
  }
}

void get_kx(real *kx, mint *kxI, const struct fourierGrid grid) {
  // Obtain the kx=(kx,ky,kz) coordinates given the multidimensional
  // kxI index. Assume the flat grid.kx array is organized as {kx,ky,kz}.
  real* kx_p = grid.kx;
  for (mint d=0; d<nDim; d++) {
    kx[d] = kx_p[kxI[d]]; 
    kx_p += grid.Nekx[d];
  }
}

// Functions that copy memory on the host.
void memcpy_real_ho(real *dest, real *src, mint numElements) {
  memcpy(dest, src, numElements*sizeof(real));
}
void memcpy_fourier_ho(fourier *dest, fourier *src, mint numElements) {
  memcpy(dest, src, numElements*sizeof(fourier));
}

// Functions that copy memory between host and device.
void memcpy_real(real *dest, real *src, mint numElements, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir != host2host)
    return memcpy_real_dev(dest, src, numElements, dir);
#endif
  memcpy_real_ho(dest, src, numElements);
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
void memcpy_fourier(void *dest, void *src, mint numElements, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir != host2host)
    return memcpy_fourier_dev(dest, src, numElements, dir);
#endif
  memcpy_fourier_ho(dest, src, numElements);
}

// Functions that copy real(Fourier)Arrays betwen host and device.
void hodevXfer_realArray(struct realArray *arr, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    memcpy_real_dev(arr->dev, arr->ho, arr->nelem, dir);
  else if (dir == device2host)
    memcpy_real_dev(arr->ho, arr->dev, arr->nelem, dir);
#endif
}
void hodevXfer_fourierArray(struct fourierArray *arr, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  if (dir == host2device)
    memcpy_fourier_dev(arr->dev, arr->ho, arr->nelem, dir);
  else if (dir == device2host)
    memcpy_fourier_dev(arr->ho, arr->dev, arr->nelem, dir);
#endif
}
