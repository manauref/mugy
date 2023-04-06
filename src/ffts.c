/* mugy: ffts

   Module with FFT operators. We do host FFTs with FFTW
   and device FFTs with cuFFT. We also anticipate 3 modes
   depending on MPI decomposition (recall data is organized
   as (z,x,y) in row-major order):
     FFTW:
       - serial or MPI decomposed along x.
       - MPI decomposed along y.
       - MPI decomposed alog x and y.
     cuFFT (don't yet know the data layout):
       - single GPU 
       - single node
       - multi node
*/

#include "mh_mpi_tools.h"
#include "mh_ffts.h"
#include "mh_data.h"
#include <string.h>   // e.g. for memcpy.

// Objects for xy FFT of a single field.
mugy_fftw_fourier *fk_xyfft_buf;  // Fourier-space buffer.
real *f_xyfft_buf;           // Real-space buffer.
mugy_fftw_plan fft_plan_xy_c2r, fft_plan_xy_r2c;
real fft_norm_xy;  // Nonunitary normalization applied in r2c transform.

void init_ffts(struct grid gridG, struct grid gridL) {

  mugy_fftw_mpi_init();

  // ....... Setup for 2D FFTs of a single field ......... //
  // Get local data size.
  const mint fftDim = 2;
  ptrdiff_t fftSizek[fftDim], fftNum, blockSizek0;
  fftSizek[0] = gridG.fG.Nekx[0];
  fftSizek[1] = gridG.fG.Nekx[1];
  fftNum      = gridL.fG.Nekx[2];
  blockSizek0 = gridL.fG.Nekx[0];
  ptrdiff_t alloc_local, local_Nekx0, local_kx0_start;
  alloc_local = mugy_fftw_mpi_local_size_many(fftDim, fftSizek, fftNum, blockSizek0, *xyComm, &local_Nekx0, &local_kx0_start);

  // Allocate buffers.
  fk_xyfft_buf = mugy_fftw_alloc_complex(alloc_local);
  f_xyfft_buf  = mugy_fftw_alloc_real(2*alloc_local);

  // Create plans.
  ptrdiff_t fftSize[fftDim], blockSize0in, blockSize0out;
  fftSize[0] = gridG.fG.dual.Nx[0];
  fftSize[1] = gridG.fG.dual.Nx[1];
  blockSize0in  = gridL.fG.dual.Nx[0];
  blockSize0out = blockSizek0;
  fft_plan_xy_r2c = mugy_fftw_mpi_plan_many_dft_r2c(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                    f_xyfft_buf, fk_xyfft_buf, *xyComm, FFTW_ESTIMATE);
//  fft_plan_xy_r2c = fftwf_mpi_plan_dft_r2c_2d(fftSize[0], fftSize[1],
//                                                  f_xyfft_buf, fk_xyfft_buf, *xyComm, FFTW_ESTIMATE);
  blockSize0in  = blockSizek0;
  blockSize0out = gridL.fG.dual.Nx[0];
  fft_plan_xy_c2r = mugy_fftw_mpi_plan_many_dft_c2r(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                    fk_xyfft_buf, f_xyfft_buf, *xyComm, FFTW_ESTIMATE);
//  fft_plan_xy_c2r = fftwf_mpi_plan_dft_c2r_2d(fftSize[0], fftSize[1],
//                                                  fk_xyfft_buf, f_xyfft_buf, *xyComm, FFTW_ESTIMATE);

  fft_norm_xy = 1./((real)gridG.fG.dual.NxyTot);
  // ....... End setup for 2D FFTs of a single field ......... //

}

void xyfft_c2r(struct realArray *fOut, struct fourierArray *fkIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return xyfft_c2r_dev(fOut, fkIn);
#endif
  
  // Copy data into buffer.
  memcpy_fourier(fk_xyfft_buf, fkIn->ho, fkIn->nelem, host2host);

  // Inverse FFT.
  mugy_fftw_mpi_execute_dft_c2r(fft_plan_xy_c2r, fk_xyfft_buf, f_xyfft_buf);

  // Copy data from buffer.
  memcpy_real(fOut->ho, f_xyfft_buf, fOut->nelem, host2host);
}

void xyfft_r2c(struct fourierArray *fkOut, struct realArray *fIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return xyfft_c2r_dev(fOut, fkIn);
#endif
  
  // Copy data into buffer.
  memcpy_real(f_xyfft_buf, fIn->ho, fIn->nelem, host2host);

  // Forward FFT.
  mugy_fftw_mpi_execute_dft_r2c(fft_plan_xy_r2c, f_xyfft_buf, fk_xyfft_buf);

  // Copy data from buffer.
  memcpy_fourier(fkOut->ho, fk_xyfft_buf, fkOut->nelem, host2host);

  // Apply the nonunitary normalization.
  scale_fourierArray(fkOut, fft_norm_xy, hostComp);
}

void terminate_ffts() {
  // Deallocate objects for xy FFT of a single field.
  mugy_fftw_destroy_plan(fft_plan_xy_c2r);
  mugy_fftw_destroy_plan(fft_plan_xy_r2c);
  mugy_fftw_free(fk_xyfft_buf);
  mugy_fftw_free(f_xyfft_buf);
}
