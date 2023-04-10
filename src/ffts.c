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

#include "mh_ffts.h"
#include <string.h>   // e.g. for memcpy.

void fft_init(struct mugy_ffts *ffts, struct mugy_grid gridG, struct mugy_grid gridL, struct mugy_comms comms) {

  mugy_fftw_mpi_init();

  // ....... Setup for 2D FFTs of a single field ......... //
  struct mugy_comms_sub *scomm = &comms.sub2d[0];
  // Get local data size.
  const mint fftDim = 2;
  ptrdiff_t fftSizek[fftDim], fftNum, blockSizek0;
  fftSizek[0] = gridG.fG.Nekx[0];
  fftSizek[1] = gridG.fG.Nekx[1];
  fftNum      = gridL.fG.Nekx[2];
  blockSizek0 = gridL.fG.Nekx[0];
  ptrdiff_t alloc_local, local_Nekx0, local_kx0_start;
  alloc_local = mugy_fftw_mpi_local_size_many(fftDim, fftSizek, fftNum, blockSizek0, scomm->comm, &local_Nekx0, &local_kx0_start);

  // Allocate buffers.
  ffts->xy.kbuf = mugy_fftw_alloc_complex(alloc_local);
  ffts->xy.rbuf = mugy_fftw_alloc_real(2*alloc_local);

  // Create plans.
  ptrdiff_t fftSize[fftDim], blockSize0in, blockSize0out;
  fftSize[0] = gridG.fG.dual.Nx[0];
  fftSize[1] = gridG.fG.dual.Nx[1];
  blockSize0in  = gridL.fG.dual.Nx[0];
  blockSize0out = blockSizek0;
  ffts->xy.plan_r2c = mugy_fftw_mpi_plan_many_dft_r2c(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                      ffts->xy.rbuf, ffts->xy.kbuf, scomm->comm, FFTW_ESTIMATE);
  blockSize0in  = blockSizek0;
  blockSize0out = gridL.fG.dual.Nx[0];
  ffts->xy.plan_c2r = mugy_fftw_mpi_plan_many_dft_c2r(fftDim, fftSize, fftNum, blockSize0in, blockSize0out,
                                                      ffts->xy.kbuf, ffts->xy.rbuf, scomm->comm, FFTW_ESTIMATE);

  ffts->xy.normFac = 1./((real)gridG.fG.dual.NxyTot);
  ffts->xy.forwardNorm = false;  // This FFT is only used for ICs given in real-space.
  // ....... End setup for 2D FFTs of a single field ......... //

}

void fft_xy_c2r(struct mugy_ffts *ffts, struct mugy_realArray *fOut, struct mugy_fourierArray *fkIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return xyfft_c2r_dev(fOut, fkIn);
#endif
  
  // Copy data into buffer.
  memcpy_fourier(ffts->xy.kbuf, fkIn->ho, fkIn->nelem, host2host);

  // Inverse FFT.
  mugy_fftw_mpi_execute_dft_c2r(ffts->xy.plan_c2r, ffts->xy.kbuf, ffts->xy.rbuf);

  // Copy data from buffer.
  memcpy_real(fOut->ho, ffts->xy.rbuf, fOut->nelem, host2host);

  // Apply the nonunitary normalization.
  if (!ffts->xy.forwardNorm)
    scale_realArray(fOut, ffts->xy.normFac, hostComp);
}

void fft_xy_r2c(struct mugy_ffts *ffts, struct mugy_fourierArray *fkOut, struct mugy_realArray *fIn, enum resource_comp res) {
#if USE_GPU
//  if (res == deviceComp)
//    return xyfft_c2r_dev(fOut, fkIn);
#endif
  
  // Copy data into buffer.
  memcpy_real(ffts->xy.rbuf, fIn->ho, fIn->nelem, host2host);

  // Forward FFT.
  mugy_fftw_mpi_execute_dft_r2c(ffts->xy.plan_r2c, ffts->xy.rbuf, ffts->xy.kbuf);

  // Copy data from buffer.
  memcpy_fourier(fkOut->ho, ffts->xy.kbuf, fkOut->nelem, host2host);

  // Apply the nonunitary normalization.
  if (ffts->xy.forwardNorm)
    scale_fourierArray(fkOut, ffts->xy.normFac, hostComp);
}

void fft_terminate(struct mugy_ffts *ffts) {
  // Deallocate objects for xy FFT of a single field.
  mugy_fftw_free(ffts->xy.kbuf);
  mugy_fftw_free(ffts->xy.rbuf);
  mugy_fftw_destroy_plan(ffts->xy.plan_c2r);
  mugy_fftw_destroy_plan(ffts->xy.plan_r2c);
}
