/* mugy: io_adios.c
 *
 * IO functions based on ADIOS.
 *
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

struct mugy_ioManager *mugy_io_init(struct mugy_comms comms) {
  // Allocate IO manager.
  struct mugy_ioManager *ioman = (struct mugy_ioManager *) malloc(sizeof(struct mugy_ioManager));

  // Initialize ADIOS IO.
  ioman->ctx = adios2_init_mpi(comms.world.comm);
  ad_check_handler(ad_check_handler, " ADIOS: Error initiating.");

  return ioman;
}

struct mugy_ad_file *mugy_io_create_mugy_array_file(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid *grid, enum mugy_datatype dtype) {
  // Create a file that will hold a mugy array defined on the grid.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = dtype == real_enum;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim], start[nDim], count[nDim];
  // The data is in z,x,y or kz,kx,ky order.
  const mint dimOrg[nDim] = {1,2,0};
  if (adf->isVarReal) {
    struct mugy_realGrid *gridG = &grid->global.deal.dual, *gridL = &grid->local.deal.dual;
    for (mint d=0; d<nDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nx", adios_mint, gridG->Nx, nDim);
    adios2_define_attribute_array(adf->io, "xMin", adios_real, gridG->xMin, nDim);
    adios2_define_attribute_array(adf->io, "xMax", adios_real, gridG->xMax, nDim);
    adios2_define_attribute_array(adf->io, "dx", adios_real, gridG->dx, nDim);
    // Variables
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_real, nDim, shape,
                                      start, count, adios2_constant_dims_true);
  } else {
    struct mugy_fourierGrid *gridG = &grid->global.deal, *gridL = &grid->local.deal;
    for (mint d=0; d<nDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nekx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nekx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nekx", adios_mint, gridG->Nekx, nDim);
    adios2_define_attribute_array(adf->io, "kxMin", adios_real, gridG->kxMin, nDim);
    // Variables
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_fourier, nDim, shape,
                                      start, count, adios2_constant_dims_true);
  }
  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  free(fnamebp);

  return adf;
}

struct mugy_ad_file *mugy_io_create_moments_file(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop, enum mugy_datatype dtype) {
  // Create a file storing real-space moments.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = dtype == real_enum;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y or s,kz,kz,ky order order.
  const mint dimOrg[nDim] = {2,3,1};
  shape[0] = (size_t)pop->global.numMomentsTot;
  start[0] = (size_t)pop->local.globalMomOff;
  count[0] = (size_t)pop->local.numMomentsTot;
  if (adf->isVarReal) {
    struct mugy_realGrid *gridG = &grid->global.deal.dual, *gridL = &grid->local.deal.dual;
    for (mint d=0; d<nDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nx", adios_mint, gridG->Nx, nDim);
    adios2_define_attribute_array(adf->io, "xMin", adios_real, gridG->xMin, nDim);
    adios2_define_attribute_array(adf->io, "xMax", adios_real, gridG->xMax, nDim);
    adios2_define_attribute_array(adf->io, "dx", adios_real, gridG->dx, nDim);
    // Variables
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_real, nDim+1, shape,
                                      start, count, adios2_constant_dims_true);
  } else {
    struct mugy_fourierGrid *gridG = &grid->global.deal, *gridL = &grid->local.deal;
    for (mint d=0; d<nDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nekx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nekx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nekx", adios_mint, gridG->Nekx, nDim);
    adios2_define_attribute_array(adf->io, "kxMin", adios_real, gridG->kxMin, nDim);
    // Variables
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_fourier, nDim+1, shape,
                                      start, count, adios2_constant_dims_true);
  }

  adios2_define_attribute(adf->io, "numSpecies", adios_mint, &pop->global.numSpecies);
  mint numMom = pop->global.numMomentsTot/pop->global.numSpecies;
  adios2_define_attribute(adf->io, "numMoments", adios_mint, &numMom);

  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  free(fnamebp);

  return adf;
}

struct mugy_ad_file *mugy_io_create_population_perp_file(struct mugy_ioManager *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop,
  enum mugy_datatype dtype, enum mugy_datatype gridtype, mint ncomp, mint zIdx) {
  // Create a file for a mugy_array holding ncomp quantities per species on an perpendicular plane.

  const mint perpDim = 2;
  bool isGridReal = gridtype == real_enum;

  // Exit if not part of this perp plane.
  mint globalOffz = isGridReal? grid->local.deal.dual.globalOff[2] : grid->local.deal.globalOff[2];;
  mint localNz    = isGridReal? grid->local.deal.dual.Nx[2] : grid->local.deal.Nekx[2];;
  if ((zIdx < globalOffz) || (globalOffz+localNz < zIdx)) return NULL;

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));
  adf->isVarReal = dtype == real_enum;

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[perpDim+1], start[perpDim+1], count[perpDim+1];
  // The data is in s,x,y or s,kz,ky order order.
  const mint dimOrg[2] = {1,2};
  shape[0] = (size_t)pop->global.numSpecies*ncomp;
  start[0] = (size_t)pop->local.globalSpecOff*ncomp;
  count[0] = (size_t)pop->local.numSpecies*ncomp;
  if (isGridReal) {
    struct mugy_realGrid *gridG = &grid->global.deal.dual, *gridL = &grid->local.deal.dual;
    for (mint d=0; d<perpDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nx", adios_mint, gridG->Nx, perpDim);
    adios2_define_attribute_array(adf->io, "xMin", adios_real, gridG->xMin, perpDim);
    adios2_define_attribute_array(adf->io, "xMax", adios_real, gridG->xMax, perpDim);
    adios2_define_attribute_array(adf->io, "dx", adios_real, gridG->dx, perpDim);
  } else {
    struct mugy_fourierGrid *gridG = &grid->global.deal, *gridL = &grid->local.deal;
    for (mint d=0; d<perpDim; d++) {
      shape[dimOrg[d]] = (size_t)gridG->Nekx[d];
      start[dimOrg[d]] = (size_t)gridL->globalOff[d];
      count[dimOrg[d]] = (size_t)gridL->Nekx[d];
    }
    // Attributes
    adios2_define_attribute_array(adf->io, "Nekx", adios_mint, gridG->Nekx, perpDim);
    adios2_define_attribute_array(adf->io, "kxMin", adios_real, gridG->kxMin, perpDim);
  }

  // Variables
  if (adf->isVarReal)
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_real, perpDim+1, shape,
                                      start, count, adios2_constant_dims_true);
  else
    adf->var = adios2_define_variable(adf->io, "globalVariable", adios_fourier, perpDim+1, shape,
                                      start, count, adios2_constant_dims_true);

  adios2_define_attribute(adf->io, "numSpecies", adios_mint, &pop->global.numSpecies);
  mint numMom = ncomp;
  adios2_define_attribute(adf->io, "numMoments", adios_mint, &numMom);

  ad_check_handler(adf->var, " ADIOS: Error defining variable.");

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  free(fnamebp);

  return adf;
}

void mugy_io_setup_files(struct mugy_ioManager *ioman, struct mugy_grid *grid, struct mugy_population *pop) {

  // List files we intend to open/write/close.
  // For some default files we prepend the type of data with this key:
  //   ra_ : real array.
  //   ka_ : fourier array.
  //   rm_ : real moments.
  //   km_ : fourier moments.
  // If more, nonstandard files wish to be added, just put the name in flist and
  // be sure to use the correct create/write functions.
  char *flist[] = {"ka_phik","km_momk"};

  ioman->numfiles = sizeof(flist)/sizeof(flist[0]);
  ioman->files = (struct mugy_ad_file **)calloc(ioman->numfiles, sizeof(struct mugy_ad_file*));

  char *dataKind = alloc_charArray_ho(3+1);
  for (mint i=0; i<ioman->numfiles; i++) {
    strncpy(dataKind, flist[i], 3);
    char *fname = alloc_charArray_ho(strlen(flist[i])-2+1);
    strcpy(fname, flist[i]+3);

    if (strcmp(dataKind, "ra_") == 0) {
      ioman->files[i] = mugy_io_create_mugy_array_file(ioman, fname, grid, real_enum);
    } else if (strcmp(dataKind, "ka_") == 0) {
      ioman->files[i] = mugy_io_create_mugy_array_file(ioman, fname, grid, fourier_enum);
    } else if (strcmp(dataKind, "rm_") == 0) {
      ioman->files[i] = mugy_io_create_moments_file(ioman, fname, grid, pop, real_enum);
    } else if (strcmp(dataKind, "km_") == 0) {
      ioman->files[i] = mugy_io_create_moments_file(ioman, fname, grid, pop, fourier_enum);
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

void mugy_io_write_mugy_array(struct mugy_ioManager *ioman, char* fname, struct mugy_ad_file *fin, struct mugy_array *arr) {
  // Write out real space array.
  adios2_error ioerr;
  struct mugy_ad_file *fh = fin == NULL ? get_fileHandle(ioman, fname) : fin;
  if (arr->ho) ioerr = adios2_put(fh->eng, fh->var, arr->ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting mugy array.");
}

void mugy_io_close_file(struct mugy_ad_file *fh) {
  // Close a file given its mugy file handle.
  adios2_error ioerr = adios2_close(fh->eng);
  ad_check_error(ioerr, " ADIOS: Error closing engine/file.");
}

void mugy_io_terminate(struct mugy_ioManager *ioman) {
  // Finalize ADIOS IO.

  // Close all files.
  for (mint i=0; i<ioman->numfiles; i++)
    mugy_io_close_file(ioman->files[i]);

  adios2_error ioerr = adios2_finalize(ioman->ctx);
  ad_check_error(ioerr, " ADIOS: Error finalizing.");

  free(ioman);
}
