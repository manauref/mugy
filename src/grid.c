/* mugy: grid.c
 *
 * Operations on having to do with the real and Fourier grids.
 *
 */

#include "mh_grid.h"
#include "mh_io_tools.h"
#include "mh_utilities.h"
#include "mh_alloc.h"
#include <stdlib.h>  // for malloc.

// Allocate mugy_grid object.
struct mugy_grid *mugy_grid_alloc() {
  struct mugy_grid *grid = (struct mugy_grid *) malloc(sizeof(struct mugy_grid));
  return grid;
}

void mugy_grid_init_global(struct mugy_grid *grid, mint rank) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  struct mugy_grid_ada *gridG = &grid->global;

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  arrPrintS_mint(gridG->deal.Nkx, nDim, " User requested  NkxG=("," ) distinct wavenumbers (absolute magnitude)\n", rank);
  for (mint d=0; d<nDim; d++) gridG->al.dual.Nx[d]  = closest_power_of_two(3*(gridG->deal.Nkx[d]-1));
  gridG->al.dual.NxTot  = prod_mint(gridG->al.dual.Nx,nDim);
  gridG->al.dual.NxyTot = prod_mint(gridG->al.dual.Nx,2);

  // Number of distinct aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) gridG->al.Nkx[d] = gridG->al.dual.Nx[d]/2+1;
  // Length of aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) gridG->al.Nekx[d] = gridG->al.Nkx[d];
  gridG->al.Nekx[0] += gridG->al.Nkx[0]-1;  // Add the negative kx's:
  gridG->al.NekxTot  = prod_mint(gridG->al.Nekx,nDim);
  gridG->al.NekxyTot = prod_mint(gridG->al.Nekx,2);

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) gridG->deal.Nkx[d] = 2*(gridG->al.Nkx[d]-1)/3+1;
  // Length of de-aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) gridG->deal.Nekx[d] = gridG->deal.Nkx[d];
  gridG->deal.Nekx[0] += gridG->deal.Nkx[0]-1;  // Add the negative kx's:
  gridG->deal.NekxTot  = prod_mint(gridG->deal.Nekx,nDim);
  gridG->deal.NekxyTot = prod_mint(gridG->deal.Nekx,2);

  // Number of cells in de-aliased real-space.
  for (mint d=0; d<nDim; d++) gridG->deal.dual.Nx[d]  = 2*(gridG->deal.Nkx[d]-1)+1;
  gridG->deal.dual.NxTot  = prod_mint(gridG->deal.dual.Nx,nDim);
  gridG->deal.dual.NxyTot = prod_mint(gridG->deal.dual.Nx,2);

  real Lx[nDim] = {2.0*M_PI/gridG->deal.kxMin[0], 2.0*M_PI/gridG->deal.kxMin[1], 2.0*M_PI/gridG->deal.kxMin[2]};

  // Length of dealised and aliased real-space cell.
  for (mint d=0; d<nDim; d++) {
    gridG->deal.dual.dx[d]  = Lx[d]/fmax(1.,(real)(gridG->deal.dual.Nx[d]-gridG->deal.dual.Nx[d] % 2));
    gridG->al.dual.dx[d] = Lx[d]/fmax(1.,(real)(gridG->al.dual.Nx[d]-gridG->al.dual.Nx[d] % 2));
  }

  // Global de-aliased real-space grids
  gridG->deal.dual.x = alloc_realArray_ho(sum_mint(gridG->deal.dual.Nx, nDim));
  real *dx = gridG->deal.dual.dx;
  mint *Nx = gridG->deal.dual.Nx; 
  mint xOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nx[d]; i++)
      gridG->deal.dual.x[i+xOff] = ((real)i)*dx[d]+(real)(1-Nx[d] % 2-1)*0.5*dx[d]-0.5*Lx[d];
    gridG->deal.dual.xMin[d] = gridG->deal.dual.x[0+xOff];
    gridG->deal.dual.xMax[d] = gridG->deal.dual.x[Nx[d]-1+xOff];
    xOff += Nx[d];
  }
  // Global aliased real-space grids (may not be needed).
  gridG->al.dual.x = alloc_realArray_ho(sum_mint(gridG->al.dual.Nx, 3));
  real *dxa = gridG->al.dual.dx;
  mint *Nxa = gridG->al.dual.Nx; 
  mint xaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nxa[d]; i++)
      gridG->al.dual.x[i+xaOff] = (real)(i)*dxa[d]+(real)(1-Nxa[d] % 2-1)*0.5*dxa[d]-0.5*Lx[d];
    gridG->al.dual.xMin[d] = gridG->al.dual.x[0+xaOff];
    gridG->al.dual.xMax[d] = gridG->al.dual.x[Nxa[d]-1+xaOff];
    xaOff += Nxa[d];
  }

  // Global dealiased k-space grids.
  for (mint d=0; d<nDim; d++) gridG->al.kxMin[d] = gridG->deal.kxMin[d];
  gridG->deal.kx  = alloc_realArray_ho(sum_mint(gridG->deal.Nekx, 3));
  real *kxMin = gridG->deal.kxMin;
  mint *Nkx = gridG->deal.Nkx; 
  mint kxOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkx[d]; i++)
      gridG->deal.kx[i+kxOff] = (real)(i)*kxMin[d];
    kxOff += gridG->deal.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (mint i=Nkx[0]; i<gridG->deal.Nekx[0]; i++)
    gridG->deal.kx[i] = -(real)(Nkx[0]-1-(i-Nkx[0]))*kxMin[0];

  // Global aliased k-space grids.
  gridG->al.kx = alloc_realArray_ho(sum_mint(gridG->al.Nekx, 3));
  real *kxaMin = gridG->al.kxMin;
  mint *Nkxa = gridG->al.Nkx; 
  mint kxaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkxa[d]; i++)
      gridG->al.kx[i+kxaOff] = (real)(i)*kxaMin[d];
    kxaOff += gridG->al.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (mint i=Nkxa[0]; i<gridG->al.Nekx[0]; i++)
    gridG->al.kx[i] = -(real)(Nkxa[0]-1-(i-Nkxa[0]))*kxaMin[0];

  r0printf("\n Proceeding with :\n", rank);
  arrPrintS_mint(gridG->deal.Nkx,      nDim, " Number of distinct de-aliased absolute wavenumbers: NkxG   =", "\n", rank);
  arrPrintS_mint(gridG->deal.Nekx,     nDim, " Length of de-aliased k-space arrays:                NekxG  =", "\n", rank);
  arrPrintS_mint(gridG->al.Nkx,        nDim, " Number of distinct aliased absolute wavenumbers:    NkxaG  =", "\n", rank);
  arrPrintS_mint(gridG->al.Nekx,       nDim, " Length of aliased k-space arrays:                   NekxaG =", "\n", rank);
  arrPrintS_mint(gridG->al.dual.Nx,    nDim, " Number of aliased real space cells:                 NxaG   =", "\n", rank);
  arrPrintS_mint(gridG->deal.dual.Nx,  nDim, " Number of de-aliased real space cells:              NxG    =", "\n", rank);

  arrPrintS_real(gridG->deal.kxMin,    nDim, " Minimum absolute magnitude of wavenumbers: kxMin    =", "\n", rank);
  arrPrintS_real(gridG->deal.kxMaxDyn, nDim, " Largest wavenumbers evolved:               kxMaxDyn =", "\n", rank);

}

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

void mugy_grid_free(struct mugy_grid *grid) {
  // Deallocate memory used by grids.
  free(grid->global.deal.dual.x);
  free(grid->global.deal.kx);
  free(grid->global.al.dual.x);
  free(grid->global.al.kx);

  free(grid->local.deal.dual.x);
  free(grid->local.deal.kx);
  free(grid->local.al.dual.x);
  free(grid->local.al.kx);
}

