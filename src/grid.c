/* mugy: grid.c
 *
 * Operations on having to do with the real and Fourier grids.
 *
 */

#include "mh_grid.h"

mint sub2lin_real(mint *xI, const struct mugy_realGrid grid) {
  // Given the nDim-dimensional index (subscript) xI return the linear index
  // in a real grid. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid.Nx[1],1,grid.NxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += xI[d]*strides[d];
  return lin;
}

mint sub2lin_fourier(mint *kxI, const struct mugy_fourierGrid grid) {
  // Given the nDim-dimensional index (subscript) kxI return the linear index
  // in a Fourier grid. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += kxI[d]*strides[d];
  return lin;
}

void lin2sub_real(mint *xI, mint lin, const struct mugy_realGrid grid) {
  // Given the linear index 'lin' in a real grid, return the nDim-dimensional
  // index (subscript) xI. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid.Nx[1],1,grid.NxyTot};
  for (mint d=0; d<nDim; d++) {
    xI[d] = lin/strides[d];
    lin -= xI[d]*strides[d];
  }
}

void lin2sub_fourier(mint *kxI, mint lin, const struct mugy_fourierGrid grid) {
  // Given the linear index 'lin' in a Fourier grid, return the nDim-dimensional
  // index (subscript) kxI. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid.Nekx[1],1,grid.NekxyTot};
  for (mint d=0; d<nDim; d++) {
    kxI[d] = lin/strides[d];
    lin -= kxI[d]*strides[d];
  }
}

void get_x(real *x, mint *xI, const struct mugy_realGrid grid) {
  // Obtain the x=(x,y,z) coordinates given the multidimensional
  // xI index. Assume the flat grid.x array is organized as {x,y,z}.
  real* x_p = grid.x;
  for (mint d=0; d<nDim; d++) {
    x[d] = x_p[xI[d]]; 
    x_p += grid.Nx[d];
  }
}

void get_kx(real *kx, mint *kxI, const struct mugy_fourierGrid grid) {
  // Obtain the kx=(kx,ky,kz) coordinates given the multidimensional
  // kxI index. Assume the flat grid.kx array is organized as {kx,ky,kz}.
  real* kx_p = grid.kx;
  for (mint d=0; d<nDim; d++) {
    kx[d] = kx_p[kxI[d]]; 
    kx_p += grid.Nekx[d];
  }
}
