/* mugy: mh_ffts_fftw_wrap.h
 * 
 * Wrap FFTW functions depending on precision.
 *
 */
#include <fftw3-mpi.h>
#include "mh_data.h"

#if USE_SINGLE_PRECISION
// ............ SINGLE PRECISION ............ // 

typedef fftwf_complex mugy_fftw_fourier;

typedef fftwf_plan mugy_fftw_plan;

static inline void mugy_fftw_mpi_init() {
  fftwf_mpi_init();
}

static inline ptrdiff_t mugy_fftw_mpi_local_size_many(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t block0, MPI_Comm comm, ptrdiff_t *local_n0, ptrdiff_t *local_0_start) {
  return fftwf_mpi_local_size_many(rnk, n, howmany, block0, comm, local_n0, local_0_start);
}

static inline real *mugy_fftw_alloc_real(size_t n) {
 return fftwf_alloc_real(n);
}
static inline fftwf_complex *mugy_fftw_alloc_complex(size_t n){
 return fftwf_alloc_complex(n);
}

static inline fftwf_plan mugy_fftw_mpi_plan_many_dft_r2c(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t iblock, ptrdiff_t oblock, real *in, fftwf_complex *out, MPI_Comm comm, unsigned flags) {
  return fftwf_mpi_plan_many_dft_r2c(rnk, n, howmany,
    iblock, oblock, in, out, comm, flags);
}
static inline fftwf_plan mugy_fftw_mpi_plan_many_dft_c2r(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t iblock, ptrdiff_t oblock, fftwf_complex *in, real *out, MPI_Comm comm, unsigned flags){
  return fftwf_mpi_plan_many_dft_c2r(rnk, n, howmany,
    iblock, oblock, in, out, comm, flags);
}

static inline void mugy_fftw_mpi_execute_dft_r2c(fftwf_plan p, real *in, fftwf_complex *out){
  fftwf_mpi_execute_dft_r2c(p, in, out);
}
static inline void mugy_fftw_mpi_execute_dft_c2r(fftwf_plan p, fftwf_complex *in, real *out){
  fftwf_mpi_execute_dft_c2r(p, in, out);
}

static inline void mugy_fftw_destroy_plan(fftwf_plan plan) {
  fftwf_destroy_plan(plan);
}

static inline void mugy_fftw_free(void *p) {
  fftwf_free(p);
}



// ............ END SINGLE PRECISION ............ // 
#else




// ............ DOUBLE PRECISION ............ // 

typedef fftw_complex mugy_fftw_fourier;

typedef fftw_plan mugy_fftw_plan;

static inline void mugy_fftw_mpi_init() {
  fftw_mpi_init();
}

static inline ptrdiff_t mugy_fftw_mpi_local_size_many(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t block0, MPI_Comm comm, ptrdiff_t *local_n0, ptrdiff_t *local_0_start) {
  return fftw_mpi_local_size_many(rnk, n, howmany, block0, comm, local_n0, local_0_start);
}

static inline double *mugy_fftw_alloc_real(size_t n) {
 return fftw_alloc_real(n);
}
static inline fftw_complex *mugy_fftw_alloc_complex(size_t n){
 return fftw_alloc_complex(n);
}

static inline fftw_plan mugy_fftw_mpi_plan_many_dft_r2c(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t iblock, ptrdiff_t oblock, double *in, fftw_complex *out, MPI_Comm comm, unsigned flags) {
  return fftw_mpi_plan_many_dft_r2c(rnk, n, howmany,
    iblock, oblock, in, out, comm, flags);
}
static inline fftw_plan mugy_fftw_mpi_plan_many_dft_c2r(int rnk, const ptrdiff_t *n, ptrdiff_t howmany,
  ptrdiff_t iblock, ptrdiff_t oblock, fftw_complex *in, double *out, MPI_Comm comm, unsigned flags){
  return fftw_mpi_plan_many_dft_c2r(rnk, n, howmany,
    iblock, oblock, in, out, comm, flags);
}

static inline void mugy_fftw_mpi_execute_dft_r2c(fftw_plan p, double *in, fftw_complex *out){
  fftw_mpi_execute_dft_r2c(p, in, out);
}
static inline void mugy_fftw_mpi_execute_dft_c2r(fftw_plan p, fftw_complex *in, double *out){
  fftw_mpi_execute_dft_c2r(p, in, out);
}

static inline void mugy_fftw_destroy_plan(fftw_plan plan) {
  fftw_destroy_plan(plan);
}

static inline void mugy_fftw_free(void *p) {
  fftw_free(p);
}

// ............ END DOUBLE PRECISION ............ // 
#endif
