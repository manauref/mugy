/* mugy: io_adios.c
   
   IO functions based on ADIOS.
*/


#include "io_adios.h"

adios2_adios *ad_ctx;  // ADIOS context used throughout our IO.

// Variable engine for outputting momk.
adios2_variable *ad_momk_var;
adios2_engine *ad_momk_eng;

void ad_check_handler(const void *handler, const char *message) {
  // Check an ADIOS handler for error.
  if (handler == NULL) abortSimulation(message);
}
void ad_check_error(const int error, const char *message) {
  // Check an ADIOS error type.
  if (error) abortSimulation(message);
}

void init_io() {
  // Initialize ADIOS IO.
  ad_ctx = adios2_init_mpi(MPI_COMM_WORLD);
  ad_check_handler(ad_ctx, " ADIOS: Error initiating.");
}

void setup_files(struct grid globalGrid, struct grid localGrid, struct population globalPop, struct population localPop) {
  // Create files for IO.

  adios2_io *ad_io = adios2_declare_io(ad_ctx, "momkWrite");
  ad_check_handler(ad_io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y order.
  shape[0] = (size_t)globalPop.numMomentsTot;
  start[0] = (size_t)localPop.globalMomOff;
  count[0] = (size_t)localPop.numMomentsTot;
  mint dimOrg[nDim] = {2,3,1};
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.Nekx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.Nekx[d];
  }

  ad_momk_var = adios2_define_variable(ad_io, "momk", adios_fourier, nDim+1, shape,
                                       start, count, adios2_constant_dims_true);
  ad_check_handler(ad_momk_var, " ADIOS: Error defining variable.");

  ad_momk_eng = adios2_open(ad_io, "momk.bp", adios2_mode_write);
  ad_check_handler(ad_momk_eng, " ADIOS: Error creating engine/opening file.");
}

void writeMoments_fourier(struct fourierMoments momkIn) {
  // Write out the moments in fourier space.
  adios2_error ioerr;
  if (momkIn.ho) ioerr = adios2_put(ad_momk_eng, ad_momk_var, momkIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier moments.");
}

void terminate_io() {
  // Finalize ADIOS IO
  adios2_error ioerr;

  ioerr = adios2_close(ad_momk_eng);
  ad_check_error(ioerr, " ADIOS: Error closing engine/file.");

  ioerr = adios2_finalize(ad_ctx);
  ad_check_error(ioerr, " ADIOS: Error finalizing.");
}
