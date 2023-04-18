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

  // Allocate space for local and global grids.
  grid->local  = (struct mugy_grid_chart *) malloc(sizeof(struct mugy_grid_chart));
  grid->global = (struct mugy_grid_chart *) malloc(sizeof(struct mugy_grid_chart));

  // Allocate for real/fourier dealiased and aliased grids.
  for (mint i=0; i<2; i++) {
    struct mugy_grid_chart *chart = i==0? grid->local : grid->global;
    chart->real      = (struct mugy_grid_basic *) malloc(sizeof(struct mugy_grid_basic));
    chart->fourier   = (struct mugy_grid_basic *) malloc(sizeof(struct mugy_grid_basic));
    chart->realAl    = (struct mugy_grid_basic *) malloc(sizeof(struct mugy_grid_basic));
    chart->fourierAl = (struct mugy_grid_basic *) malloc(sizeof(struct mugy_grid_basic));
  }

  return grid;
}

void mugy_grid_init_global(struct mugy_grid *grid, mint rank) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  struct mugy_grid_chart *gridG = grid->global;

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  arrPrintS_mint(gridG->fourier->NxNonNeg, nDim, " User requested  NkxG=("," ) distinct wavenumbers (absolute magnitude)\n", rank);
  for (mint d=0; d<nDim; d++) gridG->realAl->Nx[d]  = closest_power_of_two(3*(gridG->fourier->NxNonNeg[d]-1));
  gridG->realAl->NxTot  = prod_mint(gridG->realAl->Nx,nDim);
  gridG->realAl->NxyTot = prod_mint(gridG->realAl->Nx,2);

  // Number of distinct aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) gridG->fourierAl->NxNonNeg[d] = gridG->realAl->Nx[d]/2+1;
  // Length of aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) gridG->fourierAl->Nx[d] = gridG->fourierAl->NxNonNeg[d];
  gridG->fourierAl->Nx[0] += gridG->realAl->Nx[0] - gridG->fourierAl->NxNonNeg[0];  // Add the negative kx's.
  gridG->fourierAl->NxTot  = prod_mint(gridG->fourierAl->Nx,nDim);
  gridG->fourierAl->NxyTot = prod_mint(gridG->fourierAl->Nx,2);

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) gridG->fourier->NxNonNeg[d] = 2*(gridG->fourierAl->NxNonNeg[d]-1)/3+1;
  // Length of de-aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) gridG->fourier->Nx[d] = gridG->fourier->NxNonNeg[d];
  gridG->fourier->Nx[0] += gridG->fourier->NxNonNeg[0]-1;  // Add the negative kx's:
  gridG->fourier->NxTot  = prod_mint(gridG->fourier->Nx,nDim);
  gridG->fourier->NxyTot = prod_mint(gridG->fourier->Nx,2);

  // Number of cells in de-aliased real-space.
  for (mint d=0; d<nDim; d++) gridG->real->Nx[d]  = 2*(gridG->fourier->NxNonNeg[d]-1)+1;
  gridG->real->NxTot  = prod_mint(gridG->real->Nx,nDim);
  gridG->real->NxyTot = prod_mint(gridG->real->Nx,2);

  // Number of non-negative coordinates in real space (may not be used).
  for (mint d=0; d<nDim; d++) {
    gridG->real->NxNonNeg[d]   = mugy_max(1, gridG->real->Nx[d] / 2);    // Arbitrary, but may not be used.
    gridG->realAl->NxNonNeg[d] = mugy_max(1, gridG->realAl->Nx[d] / 2);  // Arbitrary, but may not be used.
  }

  real Lx[nDim] = {2.0*M_PI/gridG->fourier->dx[0], 2.0*M_PI/gridG->fourier->dx[1], 2.0*M_PI/gridG->fourier->dx[2]};

  // Length of dealised and aliased real-space cell.
  for (mint d=0; d<nDim; d++) {
    gridG->real->dx[d]   = Lx[d]/fmax(1.,(real)(gridG->real->Nx[d]-gridG->real->Nx[d] % 2));
    gridG->realAl->dx[d] = Lx[d]/fmax(1.,(real)(gridG->realAl->Nx[d]-gridG->realAl->Nx[d] % 2));
  }

  // Global de-aliased real-space grids
  gridG->real->x = mugy_alloc_real_ho(sum_mint(gridG->real->Nx, nDim));
  real *dx = gridG->real->dx;
  mint *Nx = gridG->real->Nx; 
  mint xOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nx[d]; i++)
      gridG->real->x[i+xOff] = ((real)i)*dx[d]+(real)(1-Nx[d] % 2-1)*0.5*dx[d]-0.5*Lx[d];
    gridG->real->xMin[d] = gridG->real->x[0+xOff];
    gridG->real->xMax[d] = gridG->real->x[Nx[d]-1+xOff];
    xOff += Nx[d];
  }
  // Global aliased real-space grids (may not be needed).
  gridG->realAl->x = mugy_alloc_real_ho(sum_mint(gridG->realAl->Nx, 3));
  real *dxa = gridG->realAl->dx;
  mint *Nxa = gridG->realAl->Nx; 
  mint xaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nxa[d]; i++)
      gridG->realAl->x[i+xaOff] = (real)(i)*dxa[d]+(real)(1-Nxa[d] % 2-1)*0.5*dxa[d]-0.5*Lx[d];
    gridG->realAl->xMin[d] = gridG->realAl->x[0+xaOff];
    gridG->realAl->xMax[d] = gridG->realAl->x[Nxa[d]-1+xaOff];
    xaOff += Nxa[d];
  }

  // Global dealiased k-space grids.
  for (mint d=0; d<nDim; d++) gridG->fourier->xMin[d] = 0.;
  gridG->fourier->x = mugy_alloc_real_ho(sum_mint(gridG->fourier->Nx, 3));
  real *dkx = gridG->fourier->dx;
  mint *Nkx = gridG->fourier->NxNonNeg; 
  mint kxOff = 0;
  // Positive Fourier modes.
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkx[d]; i++)
      gridG->fourier->x[i+kxOff] = ((real) i)*dkx[d];
    gridG->fourier->xMax[d] = gridG->fourier->x[Nkx[d]-1+kxOff];  // Arbitrary, but may not be used.
    kxOff += gridG->fourier->Nx[d];
  }
  // Negative Fourier modes in increasing order.
  for (mint i=Nkx[0]; i<gridG->fourier->Nx[0]; i++)
    gridG->fourier->x[i] = -(real)(Nkx[0]-1-(i-Nkx[0]))*dkx[0];

  // Global aliased k-space grids.
  gridG->fourierAl->x = mugy_alloc_real_ho(sum_mint(gridG->fourierAl->Nx, 3));
  real *dkxa = gridG->fourierAl->dx;
  mint *Nkxa = gridG->fourierAl->NxNonNeg; 
  mint kxaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkxa[d]; i++)
      gridG->fourierAl->x[i+kxaOff] = (real)(i)*dkxa[d];
    gridG->fourierAl->xMax[d] = gridG->fourierAl->x[Nkxa[d]-1+kxaOff];  // Arbitrary, but may not be used.
    kxaOff += gridG->fourierAl->Nx[d];
  }
  // Negative kx modes in increasing order.
  for (mint i=Nkxa[0]; i<gridG->fourierAl->Nx[0]; i++)
    gridG->fourierAl->x[i] = -(real)(Nkxa[0]-1-(i-Nkxa[0]))*dkxa[0];

  r0printf("\n Proceeding with :\n", rank);
  arrPrintS_mint(gridG->fourier->NxNonNeg,   nDim, " Number of distinct de-aliased absolute wavenumbers: NkxG   =", "\n", rank);
  arrPrintS_mint(gridG->fourier->Nx,         nDim, " Length of de-aliased k-space arrays:                NekxG  =", "\n", rank);
  arrPrintS_mint(gridG->fourierAl->NxNonNeg, nDim, " Number of distinct aliased absolute wavenumbers:    NkxaG  =", "\n", rank);
  arrPrintS_mint(gridG->fourierAl->Nx,       nDim, " Length of aliased k-space arrays:                   NekxaG =", "\n", rank);
  arrPrintS_mint(gridG->realAl->Nx,          nDim, " Number of aliased real space cells:                 NxaG   =", "\n", rank);
  arrPrintS_mint(gridG->real->Nx,            nDim, " Number of de-aliased real space cells:              NxG    =", "\n", rank);
  arrPrintS_real(gridG->fourier->dx,         nDim, " Minimum absolute magnitude of wavenumbers: kxMin    =", "\n", rank);

  // Also convenient to keep dealiased kperpSq in memory:
  gridG->fourier->xperpSq = mugy_alloc_real_ho(gridG->fourier->NxyTot);
  for (mint i=0; i<gridG->fourier->Nx[0]; i++) {
    for (mint j=0; j<gridG->fourier->Nx[1]; j++) {
      double kx = gridG->fourier->x[i];
      double ky = gridG->fourier->x[gridG->fourier->Nx[0]+j];
      gridG->fourier->xperpSq[i*gridG->fourier->Nx[1]+j] = kx*kx + ky*ky;
    }
  }

  // Set the type of each grid.
  gridG->real->type      = MUGY_REAL_GRID;
  gridG->fourier->type   = MUGY_FOURIER_GRID;
  gridG->realAl->type    = MUGY_REAL_GRID;
  gridG->fourierAl->type = MUGY_FOURIER_GRID;

}

mint mugy_grid_sub2lin_real(mint *xI, const struct mugy_grid_basic *grid) {
  // Given the nDim-dimensional index (subscript) xI return the linear index
  // in a real grid. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid->Nx[1],1,grid->NxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += xI[d]*strides[d];
  return lin;
}

mint mugy_grid_sub2lin_fourier(mint *kxI, const struct mugy_grid_basic *grid) {
  // Given the nDim-dimensional index (subscript) kxI return the linear index
  // in a Fourier grid. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid->Nx[1],1,grid->NxyTot};
  mint lin;
  for (mint d=0; d<nDim; d++) lin += kxI[d]*strides[d];
  return lin;
}

mint mugy_grid_sub2lin_perp_fourier(mint *kxI, const struct mugy_grid_basic *grid) {
  // Given the 2D perp index (subscript) kxI return the linear index
  // in a perp (kx-ky) Fourier grid. Assume row major order of (kx,ky) dimensions.
  mint strides[2] = {grid->Nx[1],1};
  mint lin;
  for (mint d=0; d<2; d++) lin += kxI[d]*strides[d];
  return lin;
}

void mugy_grid_lin2sub_real(mint *xI, mint lin, const struct mugy_grid_basic *grid) {
  // Given the linear index 'lin' in a real grid, return the nDim-dimensional
  // index (subscript) xI. We assume row major order for the (z,x,y) dimensions.
  mint strides[nDim] = {grid->Nx[1],1,grid->NxyTot};
  for (mint d=0; d<nDim; d++) {
    xI[d] = lin/strides[d];
    lin -= xI[d]*strides[d];
  }
}

void mugy_grid_lin2sub_fourier(mint *kxI, mint lin, const struct mugy_grid_basic *grid) {
  // Given the linear index 'lin' in a Fourier grid, return the nDim-dimensional
  // index (subscript) kxI. We assume row major order for the (kz,kx,ky) dimensions.
  mint strides[nDim] = {grid->Nx[1],1,grid->NxyTot};
  for (mint d=0; d<nDim; d++) {
    kxI[d] = lin/strides[d];
    lin -= kxI[d]*strides[d];
  }
}

void mugy_grid_lin2sub_fourier_perp(mint *kxI, mint lin, const struct mugy_grid_basic *grid) {
  // Given the linear index 'lin' in a perpendicular (kx-ky) Fourier grid, return the
  // 2-dimensional index (subscript) kxI. We assume row major order for the (kx,ky) dimensions.
  mint strides[2] = {grid->Nx[1],1};
  for (mint d=0; d<2; d++) {
    kxI[d] = lin/strides[d];
    lin -= kxI[d]*strides[d];
  }
}

void mugy_grid_get_x(real *x, mint *xI, const struct mugy_grid_basic *grid) {
  // Obtain the x=(x,y,z) coordinates given the multidimensional
  // xI index. Assume the flat grid.x array is organized as {x,y,z}.
  real* x_p = grid->x;
  for (mint d=0; d<nDim; d++) {
    x[d] = x_p[xI[d]]; 
    x_p += grid->Nx[d];
  }
}

void mugy_grid_get_kx(real *kx, mint *kxI, const struct mugy_grid_basic *grid) {
  // Obtain the kx=(kx,ky,kz) coordinates given the multidimensional
  // kxI index. Assume the flat grid.kx array is organized as {kx,ky,kz}.
  real* kx_p = grid->x;
  for (mint d=0; d<nDim; d++) {
    kx[d] = kx_p[kxI[d]]; 
    kx_p += grid->Nx[d];
  }
}

void mugy_grid_get_kx_perp(real *kx, mint *kxI, const struct mugy_grid_basic *grid) {
  // Obtain the kx=(kx,ky) coordinates given the multidimensional
  // kxI index. Assume the flat grid->x array is organized as {kx,ky,kz}.
  real* kx_p = grid->x;
  for (mint d=0; d<2; d++) {
    kx[d] = kx_p[kxI[d]]; 
    kx_p += grid->Nx[d];
  }
}

void mugy_grid_free(struct mugy_grid *grid) {
  // Deallocate memory used by grids.
  free(grid->global->real->x);
  free(grid->global->fourier->x);
  free(grid->global->realAl->x);
  free(grid->global->fourierAl->x);

  free(grid->local->real->x);
  free(grid->local->fourier->x);
  free(grid->local->realAl->x);
  free(grid->local->fourierAl->x);

  free(grid->global->fourier->xperpSq);
  free(grid->local->fourier->xperpSq);

  // Allocate for real/fourier dealiased and aliased grids.
  for (mint i=0; i<2; i++) {
    struct mugy_grid_chart *chart = i==0? grid->local : grid->global;
    free(chart->real     );
    free(chart->fourier  );
    free(chart->realAl   );
    free(chart->fourierAl);
  }

  free(grid->local );
  free(grid->global);

  free(grid);

}

