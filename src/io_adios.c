/* mugy: io_adios.c
   
   IO functions based on ADIOS.
*/


#include "mh_io_adios.h"

adios2_adios *ad_ctx;  // ADIOS context used throughout our IO.

// Variable engine for outputting momk.
adios2_variable *ad_momk_var;
adios2_engine *ad_momk_eng;

// Variable engine for outputting realArray.
adios2_variable *ad_arr_var;
adios2_engine *ad_arr_eng;
// Variable engine for outputting fourierArray.
adios2_variable *ad_arrk_var;
adios2_engine *ad_arrk_eng;

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

  // The data is in s,z,x,y order.
  mint dimOrg[nDim] = {2,3,1};

  // ......... File for Fourier moments ......... //
  adios2_io *ad_momk_io = adios2_declare_io(ad_ctx, "momkWrite");
  ad_check_handler(ad_momk_io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y order.
  shape[0] = (size_t)globalPop.numMomentsTot;
  start[0] = (size_t)localPop.globalMomOff;
  count[0] = (size_t)localPop.numMomentsTot;
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.Nekx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.Nekx[d];
  }

  // Attributes
  adios2_define_attribute(ad_momk_io, "numSpecies", adios_mint, &globalPop.numSpecies);
  mint numMom = globalPop.numMomentsTot/globalPop.numSpecies;
  adios2_define_attribute(ad_momk_io, "numMoments", adios_mint, &numMom);
  adios2_define_attribute_array(ad_momk_io, "Nekx", adios_mint, globalGrid.fG.Nekx, nDim);
  adios2_define_attribute_array(ad_momk_io, "kxMin", adios_real, globalGrid.fG.kxMin, nDim);
  // Variables
  ad_momk_var = adios2_define_variable(ad_momk_io, "momk", adios_fourier, nDim+1, shape,
                                       start, count, adios2_constant_dims_true);
  ad_check_handler(ad_momk_var, " ADIOS: Error defining variable.");

  ad_momk_eng = adios2_open(ad_momk_io, "momk.bp", adios2_mode_write);
  ad_check_handler(ad_momk_eng, " ADIOS: Error creating engine/opening file.");

  // ......... File for real array ......... //
  adios2_io *ad_arr_io = adios2_declare_io(ad_ctx, "arrWrite");
  ad_check_handler(ad_arr_io, " ADIOS: Error creating ADIOS r.");

  size_t shape_arr[nDim], start_arr[nDim], count_arr[nDim];
  // The data is in z,x,y order.
  for (mint d=0; d<nDim; d++)  dimOrg[d] -= 1;
  for (mint d=0; d<nDim; d++) {
    shape_arr[dimOrg[d]] = (size_t)globalGrid.fG.dual.Nx[d];
    start_arr[dimOrg[d]] = (size_t)localGrid.fG.dual.globalOff[d];
    count_arr[dimOrg[d]] = (size_t)localGrid.fG.dual.Nx[d];
  }

  // Attributes
  adios2_define_attribute_array(ad_arr_io, "Nx", adios_mint, globalGrid.fG.dual.Nx, nDim);
  adios2_define_attribute_array(ad_arr_io, "xMin", adios_real, globalGrid.fG.dual.xMin, nDim);
  adios2_define_attribute_array(ad_arr_io, "xMax", adios_real, globalGrid.fG.dual.xMax, nDim);
  adios2_define_attribute_array(ad_arr_io, "dx", adios_real, globalGrid.fG.dual.dx, nDim);
  ad_arr_var = adios2_define_variable(ad_arr_io, "arr", adios_real, nDim, shape_arr,
                                      start_arr, count_arr, adios2_constant_dims_true);
  ad_check_handler(ad_arr_var, " ADIOS: Error defining variable r.");

  ad_arr_eng = adios2_open(ad_arr_io, "arr.bp", adios2_mode_write);
  ad_check_handler(ad_arr_eng, " ADIOS: Error creating engine/opening file r.");

  // ......... File for Fourier array ......... //
  adios2_io *ad_arrk_io = adios2_declare_io(ad_ctx, "arrkWrite");
  ad_check_handler(ad_arrk_io, " ADIOS: Error creating ADIOS k.");

  size_t shape_arrk[nDim], start_arrk[nDim], count_arrk[nDim];
  // The data is in z,x,y order.
  for (mint d=0; d<nDim; d++) {
    shape_arrk[dimOrg[d]] = globalGrid.fG.Nekx[d];
    start_arrk[dimOrg[d]] = localGrid.fG.globalOff[d];
    count_arrk[dimOrg[d]] = localGrid.fG.Nekx[d];
  }

  // Attributes
  adios2_define_attribute_array(ad_arrk_io, "Nekx", adios_mint, globalGrid.fG.Nekx, nDim);
  adios2_define_attribute_array(ad_arrk_io, "kxMin", adios_real, globalGrid.fG.kxMin, nDim);
  ad_arrk_var = adios2_define_variable(ad_arrk_io, "arrk", adios_fourier, nDim, shape_arrk,
                                       start_arrk, count_arrk, adios2_constant_dims_true);
  ad_check_handler(ad_arrk_var, " ADIOS: Error defining variable k.");

  ad_arrk_eng = adios2_open(ad_arrk_io, "arrk.bp", adios2_mode_write);
  ad_check_handler(ad_arrk_eng, " ADIOS: Error creating engine/opening file k.");
}

void write_fourierMoments(struct fourierArray arrkIn) {
  // Write out Fourier space array.
  adios2_error ioerr;
  if (arrkIn.ho) ioerr = adios2_put(ad_momk_eng, ad_momk_var, arrkIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier array.");
}

void write_fourierArray(struct fourierArray arrkIn) {
  // Write out Fourier space array.
  adios2_error ioerr;
  if (arrkIn.ho) ioerr = adios2_put(ad_arrk_eng, ad_arrk_var, arrkIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier array.");
}

void write_realArray(struct realArray arrIn) {
  // Write out real space array.
  adios2_error ioerr;
  if (arrIn.ho) ioerr = adios2_put(ad_arr_eng, ad_arr_var, arrIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier array.");
}

void terminate_io() {
  // Finalize ADIOS IO
  adios2_error ioerr;

  ioerr = adios2_close(ad_momk_eng);
  ioerr = adios2_close(ad_arrk_eng);
//  ioerr = adios2_close(ad_arr_eng);
  ad_check_error(ioerr, " ADIOS: Error closing engine/file.");

  ioerr = adios2_finalize(ad_ctx);
  ad_check_error(ioerr, " ADIOS: Error finalizing.");
}
