/* mugy: io_adios.c
 *
 * IO functions based on ADIOS.
 *
 */
#include "mh_io_adios.h"
#include "mh_alloc.h"
#include <string.h>

#define MUGY_GENERATE_STRING(STRING) #STRING,
static const char *mugy_grid_types_str[] = {
  MUGY_FOREACH_GRIDTYPE(MUGY_GENERATE_STRING)
};

void ad_check_handler(const void *handler, const char *message) {
  // Check an ADIOS handler for error.
  if (handler == NULL) abortSimulation(message);
}
void ad_check_error(const int error, const char *message) {
  // Check an ADIOS error type.
  if (error) abortSimulation(message);
}

struct mugy_io *mugy_io_init(struct mugy_comms *comms) {
  // Allocate IO manager.
  struct mugy_io *ioman = (struct mugy_io *) malloc(sizeof(struct mugy_io));

  // Initialize ADIOS IO.
  ioman->ctx = adios2_init_mpi(comms->world->comm);
  ad_check_handler(ad_check_handler, " ADIOS: Error initiating.");

  ioman->world_rank = comms->world->rank;

  return ioman;
}

struct mugy_ad_file *mugy_io_create_mugy_array_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, enum mugy_data_types dtype, enum mugy_grid_types gtype) {
  // Create a file that will hold a mugy array defined on the grid.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim], start[nDim], count[nDim];
  // The data is in z,x,y or kz,kx,ky order.
  const mint dimOrg[nDim] = {1,2,0};
  struct mugy_grid_basic *gridG = gtype == MUGY_REAL_GRID ? grid->global->real : grid->global->fourier;
  struct mugy_grid_basic *gridL = gtype == MUGY_REAL_GRID ? grid->local->real :  grid->local->fourier;
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)gridG->Nx[d];
    start[dimOrg[d]] = (size_t)gridL->globalOff[d];
    count[dimOrg[d]] = (size_t)gridL->Nx[d];
  }
  // Attributes
  adios2_define_attribute_array(adf->io, "Nx", MUGY_ADIOS_MINT, gridG->Nx, nDim);
  adios2_define_attribute_array(adf->io, "dx", MUGY_ADIOS_REAL, gridG->dx, nDim);
  adios2_define_attribute(adf->io, "grid_type", MUGY_ADIOS_STRING, mugy_grid_types_str[gridG->type]);
  adios2_define_attribute_array(adf->io, "xMin", MUGY_ADIOS_REAL, gridG->xMin, nDim);
  adios2_define_attribute_array(adf->io, "xMax", MUGY_ADIOS_REAL, gridG->xMax, nDim);
  // Variables
  mint ad_dtype = dtype == MUGY_REAL? MUGY_ADIOS_REAL : MUGY_ADIOS_FOURIER;
  adf->varG = (adios2_variable **) malloc(sizeof(adios2_variable *));
  adf->varG[0] = adios2_define_variable(adf->io, "globalVariable", ad_dtype, nDim, shape,
                                        start, count, adios2_constant_dims_true);
  ad_check_handler(adf->varG[0], " ADIOS: Error defining variable.");
  if (ioman->world_rank == ioRank) {  // Only one rank should write these single numbers.
    size_t sshape[1]={1}, sstart[1]={0}, scount[1]={1};
    adf->time = adios2_define_variable(adf->io, "time", adios2_type_double, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->time, " ADIOS: Error defining variable.");
    adf->frame = adios2_define_variable(adf->io, "frame", MUGY_ADIOS_MINT, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->frame, " ADIOS: Error defining variable.");
  }

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_append);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  mugy_free_ho(fnamebp);

  return adf;
}

struct mugy_ad_file *mugy_io_create_moments_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop, enum mugy_data_types dtype, enum mugy_grid_types gtype) {
  // Create a file storing real-space moments.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[nDim+1], start[nDim+1], count[nDim+1];
  // The data is in s,z,x,y or s,kz,kz,ky order order.
  const mint dimOrg[nDim] = {2,3,1};
  shape[0] = (size_t)pop->global->numMomentsTot;
  start[0] = (size_t)pop->local->globalMomOff;
  count[0] = (size_t)pop->local->numMomentsTot;
  struct mugy_grid_basic *gridG = gtype == MUGY_REAL_GRID ? grid->global->real : grid->global->fourier;
  struct mugy_grid_basic *gridL = gtype == MUGY_REAL_GRID ? grid->local->real :  grid->local->fourier;
  for (mint d=0; d<nDim; d++) {
    shape[dimOrg[d]] = (size_t)gridG->Nx[d];
    start[dimOrg[d]] = (size_t)gridL->globalOff[d];
    count[dimOrg[d]] = (size_t)gridL->Nx[d];
  }
  // Attributes
  adios2_define_attribute_array(adf->io, "Nx", MUGY_ADIOS_MINT, gridG->Nx, nDim);
  adios2_define_attribute_array(adf->io, "dx", MUGY_ADIOS_REAL, gridG->dx, nDim);
  adios2_define_attribute(adf->io, "grid_type", MUGY_ADIOS_STRING, mugy_grid_types_str[gridG->type]);
  adios2_define_attribute_array(adf->io, "xMin", MUGY_ADIOS_REAL, gridG->xMin, nDim);
  adios2_define_attribute_array(adf->io, "xMax", MUGY_ADIOS_REAL, gridG->xMax, nDim);
  adios2_define_attribute(adf->io, "numSpecies", MUGY_ADIOS_MINT, &pop->global->numSpecies);
  mint *numMom = alloc_mintArray_ho(pop->global->numSpecies);
  for (mint s=0; s<pop->global->numSpecies; s++) numMom[s] = pop->global->pars[s].numMoments;
  adios2_define_attribute_array(adf->io, "numMoments", MUGY_ADIOS_MINT, numMom, pop->global->numSpecies);
  mugy_free_ho(numMom);
  // Variables
  mint ad_dtype = dtype == MUGY_REAL? MUGY_ADIOS_REAL : MUGY_ADIOS_FOURIER;
  adf->varG = (adios2_variable **) malloc(sizeof(adios2_variable *));
  adf->varG[0] = adios2_define_variable(adf->io, "globalVariable", ad_dtype, nDim+1, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->varG[0], " ADIOS: Error defining variable.");
  if (ioman->world_rank == ioRank) {  // Only one rank should write these single numbers.
    size_t sshape[1]={1}, sstart[1]={0}, scount[1]={1};
    adf->time = adios2_define_variable(adf->io, "time", adios2_type_double, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->time, " ADIOS: Error defining variable.");
    adf->frame = adios2_define_variable(adf->io, "frame", MUGY_ADIOS_MINT, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->frame, " ADIOS: Error defining variable.");
  }

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_append);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  mugy_free_ho(fnamebp);

  return adf;
}

struct mugy_ad_file *mugy_io_create_population_perp_file(struct mugy_io *ioman, char* fname,
  struct mugy_grid *grid, struct mugy_population *pop,
  enum mugy_data_types dtype, enum mugy_grid_types gtype, mint ncomp, mint zIdx) {
  // Create a file for a mugy_array holding ncomp quantities per species on an perpendicular plane.

  const mint perpDim = 2;

  struct mugy_grid_basic *gridG = gtype == MUGY_REAL_GRID ? grid->global->real : grid->global->fourier;
  struct mugy_grid_basic *gridL = gtype == MUGY_REAL_GRID ? grid->local->real :  grid->local->fourier;

  // Exit if not part of this perp plane.
  mint globalOffz = gridL->globalOff[2],  localNz = gridL->Nx[2];
  if ((zIdx < globalOffz) || (globalOffz+localNz < zIdx)) return NULL;

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));

  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  size_t shape[perpDim+1], start[perpDim+1], count[perpDim+1];
  // The data is in s*ncomp,x,y or s*ncomp,kz,ky order order.
  const mint dimOrg[2] = {1,2};
  shape[0] = (size_t)pop->global->numSpecies*ncomp;
  start[0] = (size_t)pop->local->globalSpecOff*ncomp;
  count[0] = (size_t)pop->local->numSpecies*ncomp;
  for (mint d=0; d<perpDim; d++) {
    shape[dimOrg[d]] = (size_t)gridG->Nx[d];
    start[dimOrg[d]] = (size_t)gridL->globalOff[d];
    count[dimOrg[d]] = (size_t)gridL->Nx[d];
  }
  // Attributes
  adios2_define_attribute_array(adf->io, "Nx", MUGY_ADIOS_MINT, gridG->Nx, perpDim);
  adios2_define_attribute_array(adf->io, "dx", MUGY_ADIOS_REAL, gridG->dx, perpDim);
  adios2_define_attribute(adf->io, "grid_type", MUGY_ADIOS_STRING, mugy_grid_types_str[gridG->type]);
  adios2_define_attribute_array(adf->io, "xMin", MUGY_ADIOS_REAL, gridG->xMin, perpDim);
  adios2_define_attribute_array(adf->io, "xMax", MUGY_ADIOS_REAL, gridG->xMax, perpDim);
  adios2_define_attribute(adf->io, "numSpecies", MUGY_ADIOS_MINT, &pop->global->numSpecies);
  mint *numMom = alloc_mintArray_ho(pop->global->numSpecies);
  for (mint s=0; s<pop->global->numSpecies; s++) numMom[s] = pop->global->pars[s].numMoments;
  adios2_define_attribute_array(adf->io, "numMoments", MUGY_ADIOS_MINT, numMom, pop->global->numSpecies);
  mugy_free_ho(numMom);
  // Variables
  mint ad_dtype = dtype == MUGY_REAL? MUGY_ADIOS_REAL : MUGY_ADIOS_FOURIER;
  adf->varG = (adios2_variable **) malloc(sizeof(adios2_variable *));
  adf->varG[0] = adios2_define_variable(adf->io, "globalVariable", ad_dtype, perpDim+1, shape,
                                    start, count, adios2_constant_dims_true);
  ad_check_handler(adf->varG[0], " ADIOS: Error defining variable.");
  if (ioman->world_rank == ioRank) {  // Only one rank should write these single numbers.
    size_t sshape[1]={1}, sstart[1]={0}, scount[1]={1};
    adf->time = adios2_define_variable(adf->io, "time", adios2_type_double, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->time, " ADIOS: Error defining variable.");
    adf->frame = adios2_define_variable(adf->io, "frame", MUGY_ADIOS_MINT, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->frame, " ADIOS: Error defining variable.");
  }

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_append);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  mugy_free_ho(fnamebp);

  return adf;
}

struct mugy_ad_file *mugy_io_create_restart_file(struct mugy_io *ioman,
  struct mugy_grid *grid, struct mugy_population *pop) {
  // Create a file storing fourier-space moments and fields to restart the simulation.

  struct mugy_ad_file *adf = (struct mugy_ad_file *) calloc(1, sizeof(struct mugy_ad_file));

  char fname[] = "mrestart";
  adf->fname = alloc_charArray_ho(strlen(fname)+1);
  strcpy(adf->fname, fname);

  adf->io = adios2_declare_io(ioman->ctx, adf->fname);
  ad_check_handler(adf->io, " ADIOS: Error creating ADIOS.");

  struct mugy_grid_basic *gridG = grid->global->fourier,  *gridL = grid->local->fourier;

  // Attributes
  adios2_define_attribute_array(adf->io, "Nx", MUGY_ADIOS_MINT, gridG->Nx, nDim);
  adios2_define_attribute_array(adf->io, "dx", MUGY_ADIOS_REAL, gridG->dx, nDim);
  adios2_define_attribute(adf->io, "grid_type", MUGY_ADIOS_STRING, mugy_grid_types_str[gridG->type]);
  adios2_define_attribute_array(adf->io, "xMin", MUGY_ADIOS_REAL, gridG->xMin, nDim);
  adios2_define_attribute_array(adf->io, "xMax", MUGY_ADIOS_REAL, gridG->xMax, nDim);
  adios2_define_attribute(adf->io, "numSpecies", MUGY_ADIOS_MINT, &pop->global->numSpecies);
  mint *numMom = alloc_mintArray_ho(pop->global->numSpecies);
  for (mint s=0; s<pop->global->numSpecies; s++) numMom[s] = pop->global->pars[s].numMoments;
  adios2_define_attribute_array(adf->io, "numMoments", MUGY_ADIOS_MINT, numMom, pop->global->numSpecies);
  mugy_free_ho(numMom);

  adf->varG = (adios2_variable **) malloc(2*sizeof(adios2_variable *));
  mint ad_dtype = MUGY_ADIOS_FOURIER;
  // Moments:
  size_t shape_mom[nDim+1], start_mom[nDim+1], count_mom[nDim+1];
  const mint dimOrg_mom[nDim] = {2,3,1};  // The data is in s,z,x,y or s,kz,kz,ky order order.
  shape_mom[0] = (size_t)pop->global->numMomentsTot;
  start_mom[0] = (size_t)pop->local->globalMomOff;
  count_mom[0] = (size_t)pop->local->numMomentsTot;
  for (mint d=0; d<nDim; d++) {
    shape_mom[dimOrg_mom[d]] = (size_t)gridG->Nx[d];
    start_mom[dimOrg_mom[d]] = (size_t)gridL->globalOff[d];
    count_mom[dimOrg_mom[d]] = (size_t)gridL->Nx[d];
  }
  adf->varG[0] = adios2_define_variable(adf->io, "momk", ad_dtype, nDim+1, shape_mom,
                                    start_mom, count_mom, adios2_constant_dims_true);
  ad_check_handler(adf->varG[0], " ADIOS: Error defining variable.");

  // Field:
  size_t shape_fld[nDim], start_fld[nDim], count_fld[nDim];
  const mint dimOrg_fld[nDim] = {1,2,0};  // The data is in z,x,y or kz,kx,ky order.
  for (mint d=0; d<nDim; d++) {
    shape_fld[dimOrg_fld[d]] = (size_t)gridG->Nx[d];
    start_fld[dimOrg_fld[d]] = (size_t)gridL->globalOff[d];
    count_fld[dimOrg_fld[d]] = (size_t)gridL->Nx[d];
  }
  adf->varG[1] = adios2_define_variable(adf->io, "phik", ad_dtype, nDim, shape_fld,
                                        start_fld, count_fld, adios2_constant_dims_true);
  ad_check_handler(adf->varG[1], " ADIOS: Error defining variable.");

  // Time, frame, time step size:
  if (ioman->world_rank == ioRank) {  // Only one rank should write these single numbers.
    size_t sshape[1]={1}, sstart[1]={0}, scount[1]={1};
    adf->time = adios2_define_variable(adf->io, "time", adios2_type_double, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->time, " ADIOS: Error defining variable.");
    adf->frame = adios2_define_variable(adf->io, "frame", MUGY_ADIOS_MINT, 1, sshape,
                                       sstart, scount, adios2_constant_dims_true);
    ad_check_handler(adf->frame, " ADIOS: Error defining variable.");
  }

  // Append .bp to file name
  char *fnamebp = alloc_charArray_ho(strlen(adf->fname)+3+1);
  strcpy(fnamebp, adf->fname);
  strcat(fnamebp,".bp");

  adf->eng = adios2_open(adf->io, fnamebp, adios2_mode_write);
  ad_check_handler(adf->eng, " ADIOS: Error creating engine/opening file.");

  mugy_free_ho(fnamebp);

  return adf;
}


void mugy_io_setup_files(struct mugy_io *ioman, struct mugy_grid *grid, struct mugy_population *pop) {

  // List files we intend to open/write/close.
  typedef struct { char* name;  bool ismom; enum mugy_data_types dtype;  enum mugy_grid_types gtype; } fop;
  static const fop flist[] = {
    {.name="phik", .ismom=false, .dtype=MUGY_FOURIER, .gtype=MUGY_FOURIER_GRID}, 
    {.name="momk", .ismom=true , .dtype=MUGY_FOURIER, .gtype=MUGY_FOURIER_GRID}, 
  };

  ioman->numfiles = sizeof(flist)/sizeof(flist[0]);
  ioman->files = (struct mugy_ad_file **)calloc(ioman->numfiles+1, sizeof(struct mugy_ad_file*));

  for (mint i=0; i<ioman->numfiles; i++) {
    if (flist[i].ismom) {
      ioman->files[i] = mugy_io_create_moments_file(ioman, flist[i].name, grid, pop, flist[i].dtype, flist[i].gtype);
    } else {
      ioman->files[i] = mugy_io_create_mugy_array_file(ioman, flist[i].name, grid, flist[i].dtype, flist[i].gtype);
    }
  }

  // Setup restart file.
  ioman->files[ioman->numfiles] = mugy_io_create_restart_file(ioman, grid, pop);

}

struct mugy_ad_file *get_fileHandle(struct mugy_io *ioman, char* fname) {
  // Given the file name, return the handle to the file in the IO manger.
  int fIdx;
  for (mint i=0; i<ioman->numfiles+1; i++) {
    if (strcmp(fname, ioman->files[i]->fname) == 0) {
      fIdx = i;  break;
    }
  }
  return ioman->files[fIdx];
}

void mugy_io_write_mugy_array(struct mugy_io *ioman, char* fname, struct mugy_ad_file *fin,
  struct mugy_array *arr, mint frame, double time) {
  // Write out a mugy array.
  adios2_error ioerr;
  struct mugy_ad_file *fh = fin == NULL ? get_fileHandle(ioman, fname) : fin;

  ioerr = adios2_put(fh->eng, fh->varG[0], arr->ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting mugy array.");

  if (ioman->world_rank == ioRank) {
    ioman->time_buf[0]  = time;
    ioerr = adios2_put(fh->eng, fh->time, ioman->time_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting time.");
    ioman->frame_buf[0] = frame;
    ioerr = adios2_put(fh->eng, fh->frame, ioman->frame_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting frame.");
  }
}

void mugy_io_write_mugy_array_step(struct mugy_io *ioman, char* fname, struct mugy_ad_file *fin,
  struct mugy_array *arr, mint frame, double time) {
  // Write out a mugy array that evolves in time.
  adios2_error ioerr;
  struct mugy_ad_file *fh = fin == NULL ? get_fileHandle(ioman, fname) : fin;

#if USE_GPU
  mugy_array_copy(arr, arr, MUGY_DEVICE2HOST);
#endif

  adios2_step_status ad_sstatus;
  ioerr = adios2_begin_step(fh->eng, adios2_step_mode_append, -1., &ad_sstatus);
  ad_check_error(ioerr, " ADIOS: Error in starting step.");

  ioerr = adios2_put(fh->eng, fh->varG[0], arr->ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting mugy array.");
  if (ioman->world_rank == ioRank) {
    ioman->time_buf[0]  = time;
    ioerr = adios2_put(fh->eng, fh->time, ioman->time_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting time.");
    ioman->frame_buf[0] = frame;
    ioerr = adios2_put(fh->eng, fh->frame, ioman->frame_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting frame.");
  }

  ioerr = adios2_end_step(fh->eng);
  ad_check_error(ioerr, " ADIOS: Error in ending step.");
}

void mugy_io_write_frame(struct mugy_io *ioman, struct mugy_population *pop,
  struct mugy_field *field, mint stepIdx, mint frame, double time) {
  // Write out a frame of all dynamical quantities of interest.

  mugy_io_write_mugy_array_step(ioman, "momk", NULL, pop->local->momk[stepIdx], frame, time);
  mugy_io_write_mugy_array_step(ioman, "phik", NULL, field->phik, frame, time);

} 

void mugy_io_write_restart(struct mugy_io *ioman, struct mugy_population *pop,
  struct mugy_field *field, struct mugy_time *time, mint stepIdx) {
  // Write out the file used to restart a simulation.
  adios2_error ioerr;
  struct mugy_ad_file *fh = get_fileHandle(ioman, "mrestart");

#if USE_GPU
  mugy_array_copy(pop->local->momk[stepIdx], pop->local->momk[stepIdx], MUGY_DEVICE2HOST);
  mugy_array_copy(field->phik, field->phik, MUGY_DEVICE2HOST);
#endif

  ioerr = adios2_put(fh->eng, fh->varG[0], pop->local->momk[stepIdx]->ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting momk.");
  ioerr = adios2_put(fh->eng, fh->varG[1], field->phik->ho, adios2_mode_deferred);
  ad_check_error(ioerr, " ADIOS: Error in putting phik.");

  if (ioman->world_rank == ioRank) {
    ioman->time_buf[0]  = time->simTime;
    ioerr = adios2_put(fh->eng, fh->time, ioman->time_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting time.");
    ioman->frame_buf[0] = time->framesOut;
    ioerr = adios2_put(fh->eng, fh->frame, ioman->frame_buf, adios2_mode_deferred);
    ad_check_error(ioerr, " ADIOS: Error in putting frame.");
  }
}

void mugy_io_close_file(struct mugy_ad_file *fh) {
  // Close a file given its mugy file handle.
  adios2_error ioerr = adios2_close(fh->eng);
  ad_check_error(ioerr, " ADIOS: Error closing engine/file.");

  free(fh->varG);  // Free array of pointers to the global variables.
}

void mugy_io_terminate(struct mugy_io *ioman) {
  // Finalize ADIOS IO.

  // Close all files.
  for (mint i=0; i<ioman->numfiles+1; i++) {
    mugy_io_close_file(ioman->files[i]);
  }

  adios2_error ioerr = adios2_finalize(ioman->ctx);
  ad_check_error(ioerr, " ADIOS: Error finalizing.");

  free(ioman);
}
