/* mugy: io_adios.c
   
   IO functions based on ADIOS.
*/


#include "mh_io_adios.h"
#include "mh_alloc.h"
#include <string.h>

void ad_check_handler(const void *handler, const char *message) {
  // Check an ADIOS handler for error.
  if (handler == NULL) abortSimulation(message);
}
void ad_check_error(const int error, const char *message) {
  // Check an ADIOS error type.
  if (error) abortSimulation(message);
}

void init_io(struct mugy_ioManager *ioman) {
  // Initialize ADIOS IO.
  ioman->ctx = adios2_init_mpi(MPI_COMM_WORLD);
  ad_check_handler(ad_check_handler, " ADIOS: Error initiating.");
}

struct mugy_ad_file *ad_create_file_realArray(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid) {
  // Create a file that will hold a real-space array.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = true;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim], start[nDim], count[nDim];
  // The data is in z,x,y order.
  mint dimOrg[nDim] = {1,2,0};
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.dual.Nx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.dual.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.dual.Nx[d];
  }

  // Attributes
  adios2_define_attribute_array(adf->io, "Nx", adios_mint, globalGrid.fG.dual.Nx, nDim);
  adios2_define_attribute_array(adf->io, "xMin", adios_real, globalGrid.fG.dual.xMin, nDim);
  adios2_define_attribute_array(adf->io, "xMax", adios_real, globalGrid.fG.dual.xMax, nDim);
  adios2_define_attribute_array(adf->io, "dx", adios_real, globalGrid.fG.dual.dx, nDim);
  // Variables
  adf->var = adios2_define_variable(adf->io, "globalVariable", adios_real, nDim, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  adf->eng = adios2_open(adf->io, strcat(adf->fname,".bp"), adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  return adf;
}

struct mugy_ad_file *ad_create_file_fourierArray(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid) {
  // Create a file that will hold a fourier-space array.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = false;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim], start[nDim], count[nDim];
  // The data is in z,x,y order.
  mint dimOrg[nDim] = {1,2,0};
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.Nekx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.Nekx[d];
  }

  // Attributes
  adios2_define_attribute_array(adf->io, "Nekx", adios_mint, globalGrid.fG.Nekx, nDim);
  adios2_define_attribute_array(adf->io, "kxMin", adios_fourier, globalGrid.fG.kxMin, nDim);
  // Variables
  adf->var = adios2_define_variable(adf->io, "globalVariable", adios_fourier, nDim, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  adf->eng = adios2_open(adf->io, strcat(adf->fname,".bp"), adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  return adf;
}

struct mugy_ad_file *ad_create_file_realMoments(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop, struct mugy_population localPop) {
  // Create a file storing real-space moments.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = true;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y order.
  shape[0] = (size_t)globalPop.numMomentsTot;
  start[0] = (size_t)localPop.globalMomOff;
  count[0] = (size_t)localPop.numMomentsTot;
  mint dimOrg[nDim+1] = {2,3,1};
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.dual.Nx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.dual.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.dual.Nx[d];
  }

  // Attributes
  adios2_define_attribute(adf->io, "numSpecies", adios_mint, &globalPop.numSpecies);
  mint numMom = globalPop.numMomentsTot/globalPop.numSpecies;
  adios2_define_attribute(adf->io, "numMoments", adios_mint, &numMom);
  adios2_define_attribute_array(adf->io, "Nx", adios_mint, globalGrid.fG.dual.Nx, nDim);
  adios2_define_attribute_array(adf->io, "xMin", adios_real, globalGrid.fG.dual.xMin, nDim);
  adios2_define_attribute_array(adf->io, "xMax", adios_real, globalGrid.fG.dual.xMax, nDim);
  adios2_define_attribute_array(adf->io, "dx", adios_real, globalGrid.fG.dual.dx, nDim);
  // Variables
  adf->var = adios2_define_variable(adf->io, "globalVariable", adios_real, nDim, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  adf->eng = adios2_open(adf->io, strcat(adf->fname,".bp"), adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  return adf;
}

struct mugy_ad_file *ad_create_file_fourierMoments(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop, struct mugy_population localPop) {
  // Create a file storing Fourier-space moments.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = false;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y order.
  shape[0] = (size_t)globalPop.numMomentsTot;
  start[0] = (size_t)localPop.globalMomOff;
  count[0] = (size_t)localPop.numMomentsTot;
  mint dimOrg[nDim+1] = {2,3,1};
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)globalGrid.fG.Nekx[d];
    start[dimOrg[d]] = (size_t)localGrid.fG.globalOff[d];
    count[dimOrg[d]] = (size_t)localGrid.fG.Nekx[d];
  }

  // Attributes
  adios2_define_attribute(adf->io, "numSpecies", adios_mint, &globalPop.numSpecies);
  mint numMom = globalPop.numMomentsTot/globalPop.numSpecies;
  adios2_define_attribute(adf->io, "numMoments", adios_mint, &numMom);
  adios2_define_attribute_array(adf->io, "Nekx", adios_mint, globalGrid.fG.Nekx, nDim);
  adios2_define_attribute_array(adf->io, "kxMin", adios_fourier, globalGrid.fG.kxMin, nDim);
  // Variables
  adf->var = adios2_define_variable(adf->io, "globalVariable", adios_fourier, nDim, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  adf->eng = adios2_open(adf->io, strcat(adf->fname,".bp"), adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  return adf;
}

void setup_files(struct mugy_ioManager *ioman, struct mugy_grid globalGrid, struct mugy_grid localGrid,
  struct mugy_population globalPop, struct mugy_population localPop) {

  // List files we intend to open/write/close.
  // For some default files we prepend the type of data with this key:
  //   ra_ : real array.
  //   ka_ : fourier array.
  //   rm_ : real moments.
  //   km_ : fourier moments.
  // If more, nonstandard files wish to be added, just put the name in flist and
  // be sure to use the correct create/write functions.
  char *flist[] = {"ra_phi"};

  ioman->numfiles = sizeof(flist)/sizeof(flist[0]);
  ioman->files = (struct mugy_ad_file **)calloc(ioman->numfiles, sizeof(struct mugy_ad_file*));

  char *dataKind = alloc_charArray_ho(3+1);
  for (mint i=0; i<ioman->numfiles; i++) {
    strncpy(dataKind, flist[i], 3);
    char *fname = alloc_charArray_ho(strlen(flist[i])-2+1);
    strcpy(fname, flist[i]+3);

    if (strcmp(dataKind, "ra_") == 0) {
      ioman->files[i] = ad_create_file_realArray(ioman, fname, globalGrid, localGrid);
    } else if (strcmp(dataKind, "ka_") == 0) {
      ioman->files[i] = ad_create_file_fourierArray(ioman, fname, globalGrid, localGrid);
    } else if (strcmp(dataKind, "rm_") == 0) {
      ioman->files[i] = ad_create_file_realMoments(ioman, fname, globalGrid, localGrid, globalPop, localPop);
    } else if (strcmp(dataKind, "fm_") == 0) {
      ioman->files[i] = ad_create_file_fourierMoments(ioman, fname, globalGrid, localGrid, globalPop, localPop);
    }
    free(fname);
  }
  free(dataKind);

}

struct mugy_ad_file *get_fileHandle(struct mugy_ioManager *ioman, char* fname) {
  // Given the file name, return the handle to the file in the IO manger.
  int fIdx;
  for (mint i=0; i<ioman->numfiles; i++) {
    if (strcmp(fname, ioman->files[i]->fname) == 0) {
      fIdx = i;  break;
    }
  }
  return ioman->files[fIdx];
}

void write_realArray(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fhIn, struct mugy_realArray arrIn) {
  // Write out real space array.
  adios2_error ioerr;
  struct mugy_ad_file *fh = fhIn == NULL ? get_fileHandle(ioman, fname) : fhIn;
  if (arrIn.ho) ioerr = adios2_put(fh->eng, fh->var, arrIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier array.");
}

void write_fourierArray(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fhIn, struct mugy_fourierArray arrIn) {
  // Write out fourier space array.
  struct mugy_ad_file *fh = fhIn == NULL ? get_fileHandle(ioman, fname) : fhIn;
  adios2_error ioerr;
  if (arrIn.ho) ioerr = adios2_put(fh->eng, fh->var, arrIn.ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting Fourier array.");
}

void io_close_file(struct mugy_ad_file *fh) {
  // Close a file given its mugy file handle.
  adios2_error ioerr = adios2_close(fh->eng);
  ad_check_error(ioerr, " ADIOS: Error closing engine/file.");
}

void io_terminate(struct mugy_ioManager *ioman) {
  // Finalize ADIOS IO.

  // Close all files.
  for (mint i=0; i<ioman->numfiles; i++)
    io_close_file(ioman->files[i]);

  adios2_error ioerr = adios2_finalize(ioman->ctx);
  ad_check_error(ioerr, " ADIOS: Error finalizing.");
}
