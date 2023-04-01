/* mugy: data_mugy
   
   Parameters and fields used throughout mugy.
*/

#include "mh_data.h"
#include "mh_data_dev.h"

struct realMoments mom, moma;
struct fourierMoments momk, momka;

fourier* getMoment_fourier(struct fourierGrid grid, struct population pop, const mint sIdx, const mint momIdx, fourier *momkIn) {
  // Return a pointer to the momIdx-th moment of the sIdx-th species in momk.
  fourier* ptrOut = momkIn;
  mint momOff = 0;
  for (mint s=0; s<sIdx; s++) momOff += pop.spec[s].numMoments;
  return ptrOut+(momOff+momIdx)*grid.NekxTot;
}

mint sub2lin_fourier(const mint *kxI, const struct fourierGrid grid) {
  // Given the nDim-dimensional index (subscript) kxI return the linear index
  // in a Fourier grid. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += kxI[d]*strides[d];
  return lin;
}

void lin2sub_fourier(mint *kxI, mint lin, const struct fourierGrid grid) {
  // Given the linear index 'lin' in a kx grid, return the nDim-dimensional
  // index (subscript) kxI. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  for (mint d=0; d<nDim; d++) {
    kxI[d] = lin/strides[d];
    lin -= kxI[d]*strides[d];
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

void memcpy_real(real *dest, real *src, mint dofs, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  memcpy_real_dev(dest, src, dofs, dir);
#endif
}
// MF 2023/03/29: Use void* because C doesn't know cuCumplex/cufourier the type.
void memcpy_fourier(void *dest, void *src, mint dofs, enum memcpy_dir_dev dir) {
#ifdef USE_GPU
  memcpy_fourier_dev(dest, src, dofs, dir);
#endif
}
