/* mugy: initialization.c
   
   Functions used to initialize the simulation.
*/

#include <string.h>   // e.g. for strcat, strlen.
#include "mh_utilities.h"
#include "mh_alloc.h"
#include "mh_comms.h"
#include "mh_io_tools.h"
#include "mh_initialization.h"
#include "mh_initialization_dev.h"
#include "mh_data.h"
#include "mh_fourier_ho.h"
#include "mh_grid.h"

// Read an variable from input file.
void readFileVar_mint(FILE *fp, const mint numElements, mint *var) {
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%d", var);
  } else {
    for (mint i=0; i<numElements; i++) fscanf(fp, "%d", &var[i]);
  }
};
void readFileVar_real(FILE *fp, const mint numElements, real *var) {
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%"fmt_real, var);
  } else {
    for (mint i=0; i<numElements; i++) fscanf(fp, "%"fmt_real, &var[i]);
  }
};

// Read species parameter composed of numElements[s] for the s-th species.
// Skip the species prior to the s-th species, and when numElements[i] is
// greater than one allocate an array for such parameter.
void readFileSpeciesPar_mint(mint **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  for (mint s=0; s<sIdx; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*d");
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%d", &(*var)[0]);
  } else {
    *var = alloc_mintArray_ho(numElements[sIdx]);
    for (mint i=0; i<numElements[sIdx]; i++) fscanf(fp, "%d", &(*var)[i]);
  }
  // Skip species after sIdx (presumably will be read later).
  for (mint s=sIdx+1; s<numSpecies; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*d");
  }
}
void readFileSpeciesPar_real(real **var, FILE *fp, const mint sIdx, const mint numSpecies, const mint *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  for (mint s=0; s<sIdx; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*"fmt_real);
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%"fmt_real, &(*var)[0]);
  } else {
    *var = alloc_realArray_ho(numElements[sIdx]);
    for (mint i=0; i<numElements[sIdx]; i++) fscanf(fp, "%"fmt_real, &(*var)[i]);
  }
  // Skip species after sIdx (presumably will be read later).
  for (mint s=sIdx+1; s<numSpecies; s++) {
    for (mint i=0; i<numElements[s]; i++) fscanf(fp, "%*"fmt_real);
  }
}

void read_inputFile(const char *fileNameIn, struct mugy_grid *grid, struct mugy_timeSetup *time,
                    struct mugy_population *pop, struct mugy_fieldParameters *field, mint rank) {
  // Read input values from input file.

  if (rank == ioRank) {  // Only ioRank reads from input file.
    printf(" Reading inputs from %s\n\n",fileNameIn);

    FILE *file_p = fopen(fileNameIn, "r");  // Open for read only.

    fscanf(file_p, "%*s");  // &space.
    readFileVar_mint(file_p, nDim, &grid->fG.Nkx[0]);
    readFileVar_real(file_p, nDim, &grid->fG.kxMin[0]);
    readFileVar_real(file_p, nDim, grid->fG.kxMaxDyn);
    readFileVar_mint(file_p, nDim, grid->mpiProcs);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &time.
    readFileVar_real(file_p, 1, &time->dt);
    readFileVar_real(file_p, 1, &time->endTime);
    readFileVar_mint(file_p, 1, &time->nFrames);
    // ARKODE parameters:
    readFileVar_mint(file_p, 1, &time->ark_kySplit);
    readFileVar_mint(file_p, 1, &time->ark_fastTableExp);
    readFileVar_mint(file_p, 1, &time->ark_fastTableImp);
    readFileVar_mint(file_p, 1, &time->ark_slowTable);
    readFileVar_real(file_p, 1, &time->ark_dtFast);
    readFileVar_real(file_p, 1, &time->ark_rtol);
    readFileVar_real(file_p, 1, &time->ark_atol);
    readFileVar_mint(file_p, 1, &time->ark_ewtScaling);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &species.
    readFileVar_mint(file_p, 1, &pop->numSpecies);
    readFileVar_mint(file_p, 1, &pop->mpiProcs);
    pop->spec = (struct mugy_species*) calloc(pop->numSpecies, sizeof(struct mugy_species));
    mint *specNumMoms = (mint*) calloc(pop->numSpecies, sizeof(mint));
    mint *specOnes    = (mint*) calloc(pop->numSpecies, sizeof(mint));
    mint *specnDim    = (mint*) calloc(pop->numSpecies, sizeof(mint));
    readFileVar_mint(file_p, pop->numSpecies, specNumMoms);
    mint filePos = ftell(file_p);
    for (mint s=0; s<pop->numSpecies; s++) {
      pop->spec[s].numMoments = specNumMoms[s];
      specOnes[s] = 1;
      specnDim[s] = nDim;
    }
    real* real_p; real** real_pp; mint* mint_p;
    for (mint s=0; s<pop->numSpecies; s++) {
      fseek(file_p, filePos, SEEK_SET);
      real_p  =    &pop->spec[s].qCharge; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =     &pop->spec[s].muMass; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =        &pop->spec[s].tau; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =       &pop->spec[s].omSt; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =        &pop->spec[s].omd; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =      &pop->spec[s].delta; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =  &pop->spec[s].deltaPerp; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =        &pop->spec[s].eta; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_pp =      &pop->spec[s].alpha; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specNumMoms);
      real_pp =         &pop->spec[s].nu; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specNumMoms);
      real_p  =     &pop->spec[s].delta0; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_pp = &pop->spec[s].hDiffOrder; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specnDim   );
      real_pp =      &pop->spec[s].hDiff; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specnDim   );
      real_pp =   &pop->spec[s].kDiffMin; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specnDim   );
      mint_p  =       &pop->spec[s].icOp; readFileSpeciesPar_mint(&mint_p, file_p, s, pop->numSpecies, specOnes   );
      real_pp =    &pop->spec[s].initAux; readFileSpeciesPar_real(real_pp, file_p, s, pop->numSpecies, specnDim   );
      real_p  =      &pop->spec[s].initA; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
      real_p  =     &pop->spec[s].noiseA; readFileSpeciesPar_real(&real_p, file_p, s, pop->numSpecies, specOnes   );
    }
    free(specnDim); free(specOnes); free(specNumMoms);
    fscanf(file_p, "%*s");  // /.
  
    fscanf(file_p, "%*s");  // &field.
    readFileVar_real(file_p, 1, &field->lambdaD);
    readFileVar_mint(file_p, 1, &field->pade);
    readFileVar_mint(file_p, 1, &field->icOp);
    fscanf(file_p, "%*s");  // /.

    fclose(file_p);
  }

  // Broadcast to other processes.
  MPI_Bcast(&grid->fG.Nkx     , nDim, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->fG.kxMin   , nDim, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->fG.kxMaxDyn, nDim, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->mpiProcs   , nDim, mpi_mint, ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&time->dt              , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->endTime         , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->nFrames         , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_kySplit     , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableExp, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableImp, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_slowTable   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_dtFast      , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_rtol        , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_atol        , 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_ewtScaling  , 1, mpi_mint, ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&pop->numSpecies, 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&pop->mpiProcs  , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  if (rank != ioRank) pop->spec = (struct mugy_species*) calloc(pop->numSpecies, sizeof(struct mugy_species));
  for (mint s=0; s<pop->numSpecies; s++) {
    MPI_Bcast(&pop->spec[s].numMoments,                       1, mpi_mint, ioRank, MPI_COMM_WORLD);
    if (rank != ioRank) {
      pop->spec[s].alpha      = alloc_realArray_ho(pop->spec[s].numMoments);
      pop->spec[s].nu         = alloc_realArray_ho(pop->spec[s].numMoments);
      pop->spec[s].hDiffOrder = alloc_realArray_ho(nDim);
      pop->spec[s].hDiff      = alloc_realArray_ho(nDim);
      pop->spec[s].kDiffMin   = alloc_realArray_ho(nDim);
      pop->spec[s].initAux    = alloc_realArray_ho(nDim);
    }
    MPI_Bcast(&pop->spec[s].qCharge   ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].muMass    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].tau       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].omSt      ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].omd       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].delta     ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].deltaPerp ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].eta       ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].alpha      , pop->spec[s].numMoments, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].nu         , pop->spec[s].numMoments, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].delta0    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].hDiffOrder ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].hDiff      ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].kDiffMin   ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].icOp      ,                       1, mpi_mint, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].initAux    ,                    nDim, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].initA     ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].noiseA    ,                       1, mpi_real, ioRank, MPI_COMM_WORLD);
  }

  MPI_Bcast(&field->lambdaD, 1, mpi_real, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->pade   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->icOp   , 1, mpi_mint, ioRank, MPI_COMM_WORLD);

}

void read_inputs(mint argc, char *argv[], struct mugy_ioSetup *ioSet, struct mugy_grid *grid, struct mugy_timeSetup *time,
                 struct mugy_population *pop, struct mugy_fieldParameters *field, mint rank) {
  // Read inputs from command line arguments and input file.

  // Check for commandline arguments.
  // Currently we expect two arguments in this order:
  //   1) name of the input file.
  //   2) absolute address of the output directory.
  // For restarts add the restart directory as a 3rd argument.
  if (argc < 3) {
    printf("\n --> Not enough inputs. Need input file and output directory as command line arguments.\n");
    abortSimulation(" TERMINATING ");
  } else {
    ioSet->inputFile = argv[1];
    ioSet->outputDir = argv[2];
    ioSet->isRestart   = false;
    ioSet->outToOldDir = false;
    if (argc > 3) {  // Simulation is a restart of a previous one.
      ioSet->restartDir = argv[3];
      ioSet->isRestart  = true;
      // Output to the same directory as the previous run?
      char *checkFile = malloc(strlen(ioSet->outputDir)+strlen("phik.bp"));
      checkFile[0] = '\0';
      strcat(strcat(checkFile,ioSet->outputDir),"phik.bp");
      ioSet->outToOldDir = fileExists(checkFile);
      free(checkFile);
    }
  }

  read_inputFile(ioSet->inputFile, grid, time, pop, field, rank);

  // Set the total number of moments.
  pop->numMomentsTot = 0;
  for (mint s=0; s<pop->numSpecies; s++) pop->numMomentsTot += pop->spec[s].numMoments;

  // Set the moments pointer in the global population to NULL (we don't store global moments).
  pop->momk = NULL;

}

void init_global_grids(struct mugy_grid *globalGrid, mint rank) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  arrPrintS_mint(globalGrid->fG.Nkx, nDim, " User requested  NkxG=("," ) distinct wavenumbers (absolute magnitude)\n", rank);
  for (mint d=0; d<nDim; d++) globalGrid->fGa.dual.Nx[d]  = closest_power_of_two(3*(globalGrid->fG.Nkx[d]-1));
  globalGrid->fGa.dual.NxTot  = prod_mint(globalGrid->fGa.dual.Nx,nDim);
  globalGrid->fGa.dual.NxyTot = prod_mint(globalGrid->fGa.dual.Nx,2);

  // Number of distinct aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) globalGrid->fGa.Nkx[d] = globalGrid->fGa.dual.Nx[d]/2+1;
  // Length of aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) globalGrid->fGa.Nekx[d] = globalGrid->fGa.Nkx[d];
  globalGrid->fGa.Nekx[0] += globalGrid->fGa.Nkx[0]-1;  // Add the negative kx's:
  globalGrid->fGa.NekxTot  = prod_mint(globalGrid->fGa.Nekx,nDim);
  globalGrid->fGa.NekxyTot = prod_mint(globalGrid->fGa.Nekx,2);

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  for (mint d=0; d<nDim; d++) globalGrid->fG.Nkx[d] = 2*(globalGrid->fGa.Nkx[d]-1)/3+1;
  // Length of de-aliased arrays along kx and ky.
  for (mint d=0; d<nDim; d++) globalGrid->fG.Nekx[d] = globalGrid->fG.Nkx[d];
  globalGrid->fG.Nekx[0] += globalGrid->fG.Nkx[0]-1;  // Add the negative kx's:
  globalGrid->fG.NekxTot  = prod_mint(globalGrid->fG.Nekx,nDim);
  globalGrid->fG.NekxyTot = prod_mint(globalGrid->fG.Nekx,2);

  // Number of cells in de-aliased real-space.
  for (mint d=0; d<nDim; d++) globalGrid->fG.dual.Nx[d]  = 2*(globalGrid->fG.Nkx[d]-1)+1;
  globalGrid->fG.dual.NxTot  = prod_mint(globalGrid->fG.dual.Nx,nDim);
  globalGrid->fG.dual.NxyTot = prod_mint(globalGrid->fG.dual.Nx,2);

  real Lx[nDim] = {2.0*M_PI/globalGrid->fG.kxMin[0], 2.0*M_PI/globalGrid->fG.kxMin[1], 2.0*M_PI/globalGrid->fG.kxMin[2]};

  // Length of dealised and aliased real-space cell.
  for (mint d=0; d<nDim; d++) {
    globalGrid->fG.dual.dx[d]  = Lx[d]/fmax(1.,(real)(globalGrid->fG.dual.Nx[d]-globalGrid->fG.dual.Nx[d] % 2));
    globalGrid->fGa.dual.dx[d] = Lx[d]/fmax(1.,(real)(globalGrid->fGa.dual.Nx[d]-globalGrid->fGa.dual.Nx[d] % 2));
  }

  // Global de-aliased real-space grids
  globalGrid->fG.dual.x = alloc_realArray_ho(sum_mint(globalGrid->fG.dual.Nx, nDim));
  real *dx = globalGrid->fG.dual.dx;
  mint *Nx = globalGrid->fG.dual.Nx; 
  mint xOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nx[d]; i++)
      globalGrid->fG.dual.x[i+xOff] = ((real)i)*dx[d]+(real)(1-Nx[d] % 2-1)*0.5*dx[d]-0.5*Lx[d];
    globalGrid->fG.dual.xMin[d] = globalGrid->fG.dual.x[0+xOff];
    globalGrid->fG.dual.xMax[d] = globalGrid->fG.dual.x[Nx[d]-1+xOff];
    xOff += Nx[d];
  }
  // Global aliased real-space grids (may not be needed).
  globalGrid->fGa.dual.x = alloc_realArray_ho(sum_mint(globalGrid->fGa.dual.Nx, 3));
  real *dxa = globalGrid->fGa.dual.dx;
  mint *Nxa = globalGrid->fGa.dual.Nx; 
  mint xaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nxa[d]; i++)
      globalGrid->fGa.dual.x[i+xaOff] = (real)(i)*dxa[d]+(real)(1-Nxa[d] % 2-1)*0.5*dxa[d]-0.5*Lx[d];
    globalGrid->fGa.dual.xMin[d] = globalGrid->fGa.dual.x[0+xaOff];
    globalGrid->fGa.dual.xMax[d] = globalGrid->fGa.dual.x[Nxa[d]-1+xaOff];
    xaOff += Nxa[d];
  }

  // Global dealiased k-space grids.
  for (mint d=0; d<nDim; d++) globalGrid->fGa.kxMin[d] = globalGrid->fG.kxMin[d];
  globalGrid->fG.kx  = alloc_realArray_ho(sum_mint(globalGrid->fG.Nekx, 3));
  real *kxMin = globalGrid->fG.kxMin;
  mint *Nkx = globalGrid->fG.Nkx; 
  mint kxOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkx[d]; i++)
      globalGrid->fG.kx[i+kxOff] = (real)(i)*kxMin[d];
    kxOff += globalGrid->fG.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (mint i=Nkx[0]; i<globalGrid->fG.Nekx[0]; i++)
    globalGrid->fG.kx[i] = -(real)(Nkx[0]-1-(i-Nkx[0]))*kxMin[0];

  // Global aliased k-space grids.
  globalGrid->fGa.kx = alloc_realArray_ho(sum_mint(globalGrid->fGa.Nekx, 3));
  real *kxaMin = globalGrid->fGa.kxMin;
  mint *Nkxa = globalGrid->fGa.Nkx; 
  mint kxaOff = 0;
  for (mint d=0; d<nDim; d++) {
    for (mint i=0; i<Nkxa[d]; i++)
      globalGrid->fGa.kx[i+kxaOff] = (real)(i)*kxaMin[d];
    kxaOff += globalGrid->fGa.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (mint i=Nkxa[0]; i<globalGrid->fGa.Nekx[0]; i++)
    globalGrid->fGa.kx[i] = -(real)(Nkxa[0]-1-(i-Nkxa[0]))*kxaMin[0];

  r0printf("\n Proceeding with :\n", rank);
  arrPrintS_mint(globalGrid->fG.Nkx,      nDim, " Number of distinct de-aliased absolute wavenumbers: NkxG   =", "\n", rank);
  arrPrintS_mint(globalGrid->fG.Nekx,     nDim, " Length of de-aliased k-space arrays:                NekxG  =", "\n", rank);
  arrPrintS_mint(globalGrid->fGa.Nkx,     nDim, " Number of distinct aliased absolute wavenumbers:    NkxaG  =", "\n", rank);
  arrPrintS_mint(globalGrid->fGa.Nekx,    nDim, " Length of aliased k-space arrays:                   NekxaG =", "\n", rank);
  arrPrintS_mint(globalGrid->fGa.dual.Nx, nDim, " Number of aliased real space cells:                 NxaG   =", "\n", rank);
  arrPrintS_mint(globalGrid->fG.dual.Nx,  nDim, " Number of de-aliased real space cells:              NxG    =", "\n", rank);

  arrPrintS_real(globalGrid->fG.kxMin,    nDim, " Minimum absolute magnitude of wavenumbers: kxMin    =", "\n", rank);
  arrPrintS_real(globalGrid->fG.kxMaxDyn, nDim, " Largest wavenumbers evolved:               kxMaxDyn =", "\n", rank);

}

void allocate_dynfields(struct mugy_grid localGrid, struct mugy_population *localPop) {
  // Allocate various fields needed.
#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Allocate moments vector needed for time stepping.
  localPop->momk = (struct mugy_array*) calloc(TIME_STEPPER_NUM_FIELDS, sizeof(struct mugy_array));
  for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
    alloc_fourierMoments(&localPop->momk[s], localGrid.fG, *localPop, onResource);

}

void set_initialCondition(struct mugy_grid globalGrid, struct mugy_grid localGrid, struct mugy_population globalPop,
  struct mugy_population *localPop, struct mugy_ffts *fftMan, struct mugy_ioManager *ioman) {
  // Impose the initial conditions on the moments and thoe potential.

  struct mugy_array momk = localPop->momk[0]; // Put ICs in first stepper field.

  // NOTE: For now assume initialOp is the same for all species.
  mint initialOp = localPop->spec[0].icOp; 

  if (initialOp == 0) {
    // Initialize in real space and transform to Fourier.
    struct mugy_realGrid *grid = &localGrid.fG.dual;
    struct mugy_array momIC;
    alloc_realMoments(&momIC, *grid, *localPop, hostAndDeviceMem);

    for (mint s=0; s<localPop->numSpecies; s++) {
      real initA = localPop->spec[s].initA;

      real *den_p  = getMoment_real(*grid, *localPop, s, denIdx, momIC.ho);  // Get density of species s.
      real *temp_p = getMoment_real(*grid, *localPop, s, tempIdx, momIC.ho);  // Get temperature of species s.

      for (mint linIdx=0; linIdx<grid->NxTot; linIdx++) {
        mint xIdx[nDim];
        lin2sub_real(&xIdx[0], linIdx, *grid);  // Convert linear index to multidimensional x index.
        real x[nDim];
        get_x(&x[0], xIdx, *grid);

        // Initial density: a superposition of sines and cosines.
        den_p[0] = 0.;
        double kx = 0.1;
        for (int i; i<localGrid.fG.Nkx[0]; i++) {
          kx += i*0.2;
          double ky = 0.3;
          for (int j; j<localGrid.fG.Nkx[1]; j++) {
            ky += i*0.1;
            den_p[0] += initA*sin(kx*x[0])*cos(ky*x[1]);
          }
        }
        den_p++;
        
        // Initial temperature = 0.
        temp_p[0] = 0.;
        temp_p++;
      }
    }

    // Copy initialized moments from host to device.
    mugy_array_hodevXfer(&momIC, host2device);

    // FFT moments.
    //fft_moments_r2c(fftMan, &momk, &momIC, deviceComp)

    //......................................................
    // FFT test
    struct mugy_array fxy_r, fxy_k;
    mugy_array_alloc(&fxy_r, real_enum, grid->NxTot, hostMem);
    mugy_array_alloc(&fxy_k, fourier_enum, localGrid.fG.NekxTot, hostMem);

    // Assign real array to a linear combo of sines and cosines.
    real *fxy_rp = fxy_r.ho;
    for (mint linIdx=0; linIdx<grid->NxTot; linIdx++) {
      real initA = localPop->spec[0].initA;
      mint xIdx[nDim];
      lin2sub_real(&xIdx[0], linIdx, *grid);  // Convert linear index to multidimensional x index.
      real x[nDim];
      get_x(&x[0], xIdx, *grid);

      fxy_rp[0] = 0.;
      double kx = localGrid.fG.kxMin[0];
      double ky = localGrid.fG.kxMin[1];
      fxy_rp[0] += initA*sin(kx*x[0])*cos(ky*x[1]);
      fxy_rp++;
    }

    fft_xy_r2c(fftMan, &fxy_k, &fxy_r, hostComp);
    fft_xy_c2r(fftMan, &fxy_r, &fxy_k, hostComp);

    struct mugy_ad_file *fh = ad_create_mugy_array_file(ioman, "arr", globalGrid, localGrid, real_enum);
    write_mugy_array(NULL, "arr", fh, fxy_r);
    io_close_file(fh);

    mugy_array_free(&fxy_r, hostMem);
    mugy_array_free(&fxy_k, hostMem);
    //......................................................

    mugy_array_free(&momIC, hostMem);

  } else if (initialOp == 1) {
    // Initialize with a k-spce power law.
    real *kxMin = &localGrid.fG.kxMin[0];

    for (mint s=0; s<localPop->numSpecies; s++) {
      real initA    = localPop->spec[s].initA;
      real *initAux = &localPop->spec[s].initAux[0];

      fourier *den_p  = getMoment_fourier(localGrid.fG, *localPop, s, denIdx, momk.ho);  // Get density of species s.
      fourier *temp_p = getMoment_fourier(localGrid.fG, *localPop, s, tempIdx, momk.ho);  // Get temperature of species s.

      for (mint linIdx=0; linIdx<localGrid.fG.NekxTot; linIdx++) {
        mint kxIdx[nDim];
        lin2sub_fourier(&kxIdx[0], linIdx, localGrid.fG);  // Convert linear index to multidimensional kx index.
        real kx[nDim];
        get_kx(&kx[0], kxIdx, localGrid.fG);
  
        // Set density to a power-law in k-space.
        den_p[0] = initA*(pow((kxMin[0]+fabs(kx[0]))/kxMin[0],initAux[0]))
                        *(pow((kxMin[1]+fabs(kx[1]))/kxMin[1],initAux[1]));
        den_p++;
  
        // Set the initial temperature (fluctuations) to zero.
        temp_p[0] = 0.;
        temp_p++;
      };
    };

    // Copy initialized moments from host to device.
    mugy_array_hodevXfer(&momk, host2device);
  }

}

void init_all(mint argc, char *argv[], struct mugy_comms *comms, struct mugy_ioManager *ioman,
  struct mugy_grid *gridG, struct mugy_grid *gridL, struct mugy_timeSetup *timePars,
  struct mugy_population *popG, struct mugy_population *popL,
  struct mugy_fieldParameters *fieldPars, struct mugy_ffts *fftMan) {
  // Run the full initialization.

  // Read inputs (from command line arguments and input file).
  read_inputs(argc, argv, &ioman->setup, gridG, timePars, popG, fieldPars, comms->world.rank);

#ifdef USE_GPU
  // Initialize devices (GPUs) if any.
  init_dev(comms->world.rank);
#endif

  init_io(ioman);  // Initialize IO interface.

  comms_sub_init(comms, *gridG, *popG);  // Initialize sub-communicators.

  // Set the number of cells in Fourier space and aliased real space.
  init_global_grids(gridG, comms->world.rank);

  // Decompose the x,y,z,s domains amongst MPI processes.
  distributeDOFs(*comms, *gridG, *popG, gridL, popL);

  allocate_dynfields(*gridL, popL);  // Allocate dynamic fields.

  fft_init(fftMan, *gridG, *gridL, *comms);  // Initialize FFT infrastructure.

  setup_files(ioman, *gridG, *gridL, *popG, *popL);  // Setup IO files.

  set_initialCondition(*gridG, *gridL, *popG, popL, fftMan, ioman);  // Impose ICs.
}

void free_fields() {
  // Deallocate fields.
//#ifdef USE_GPU
//  enum resource_mem onResource = hostAndDeviceMem;
//#else
//  enum resource_mem onResource = hostMem;
//#endif
}

void free_grid(struct mugy_grid *grid) {
  // Deallocate memory used by grids.
  free(grid->fG.dual.x);
  free(grid->fG.kx);
  free(grid->fGa.dual.x);
  free(grid->fGa.kx);
}

void free_population(struct mugy_population *pop) {
  // Deallocate memory used in species struct.
  for (mint s=0; s<pop->numSpecies; s++) {
    free(pop->spec[s].alpha);
    free(pop->spec[s].nu);
    free(pop->spec[s].hDiffOrder);
    free(pop->spec[s].hDiff);
    free(pop->spec[s].kDiffMin);
    free(pop->spec[s].initAux);
  }
  free(pop->spec);

#ifdef USE_GPU
  enum resource_mem onResource = hostAndDeviceMem;
#else
  enum resource_mem onResource = hostMem;
#endif

  // Free moments vector.
  if (pop->momk) {
    for (mint s=0; s<TIME_STEPPER_NUM_FIELDS; s++)
      mugy_array_free(&pop->momk[s], onResource);
  }
}
