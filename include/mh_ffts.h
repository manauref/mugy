/* mugy: mh_ffts.h

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

#ifndef MUGY_FFTS
#define MUGY_FFTS

#include "mh_userFLAGS.h"
#include "mh_fftw_wrap.h"

void init_ffts(struct mugy_grid gridG, struct mugy_grid gridL);

void xyfft_c2r(struct mugy_realArray *fOut, struct mugy_fourierArray *fkIn, enum resource_comp res);
void xyfft_r2c(struct mugy_fourierArray *fkOut, struct mugy_realArray *fIn, enum resource_comp res);

void terminate_ffts();

#endif
