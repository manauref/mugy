/* mugy
   
   Functions used to initialize the simulation.
*/
#include "initialization.h"

void readFileVar_int(FILE *fp, const int numElements, int *var) {
  // Read an int variable from input file.
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%d", var);
  } else {
    for (int i=0; i<numElements; i++) fscanf(fp, "%d", &var[i]);
  }
};
void readFileVar_real(FILE *fp, const int numElements, real *var) {
  // Read a real variable from input file.
  fscanf(fp, "%*s = ");
  if (numElements == 1) {
    fscanf(fp, "%"SCNfREAL, var);
  } else {
    for (int i=0; i<numElements; i++) fscanf(fp, "%"SCNfREAL, &var[i]);
  }
};

// Read species parameter composed of numElements[s] for the s-th species.
// Skip the species prior to the s-th species, and when numElements[i] is
// greater than one allocate an array for such parameter.
void readFileSpeciesPar_int(int **var, FILE *fp, const int sIdx, const int *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  for (int s=0; s<sIdx; s++) {
    for (int i=0; i<numElements[s]; i++) fscanf(fp, "%*d");
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%d", &(*var)[0]);
  } else {
    *var = alloc_intArray(numElements[sIdx]);
    for (int i=0; i<numElements[sIdx]; i++) fscanf(fp, "%d", &(*var)[i]);
  }
}
void readFileSpeciesPar_real(real **var, FILE *fp, const int sIdx, const int *numElements) {
  fscanf(fp, "%*s = ");
  // Skip species prior to sIdx (presumably already read).
  for (int s=0; s<sIdx; s++) {
    for (int i=0; i<numElements[s]; i++) fscanf(fp, "%*"SCNfREAL);
  }
  // Read this species' parameter.
  if (numElements[sIdx] == 1) {
    fscanf(fp, "%"SCNfREAL, &(*var)[0]);
  } else {
    *var = alloc_realArray(numElements[sIdx]);
    for (int i=0; i<numElements[sIdx]; i++) fscanf(fp, "%"SCNfREAL, &(*var)[i]);
  }
}

void read_inputFile(const char *fileNameIn, struct grid *grid, struct timeSetup *time,
                    struct population *pop, struct fieldParameters *field) {
  // Read input values from input file.

  if (myRank == ioRank) {  // Only ioRank reads from input file.
    printf(" Reading inputs from %s\n\n",fileNameIn);

    FILE *file_p = fopen(fileNameIn, "r");  // Open for read only.

    fscanf(file_p, "%*s");  // &space.
    readFileVar_int( file_p, nDim, &grid->fG.Nkx[0]);
    readFileVar_real(file_p, nDim, &grid->fG.kxMin[0]);
    readFileVar_real(file_p, nDim, grid->fG.kxMaxDyn);
    readFileVar_int( file_p, nDim, grid->mpiProcs);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &time.
    readFileVar_real(file_p, 1, &time->dt);
    readFileVar_real(file_p, 1, &time->endTime);
    readFileVar_int( file_p, 1, &time->nFrames);
    // ARKODE parameters:
    readFileVar_int( file_p, 1, &time->ark_kySplit);
    readFileVar_int( file_p, 1, &time->ark_fastTableExp);
    readFileVar_int( file_p, 1, &time->ark_fastTableImp);
    readFileVar_int( file_p, 1, &time->ark_slowTable);
    readFileVar_real(file_p, 1, &time->ark_dtFast);
    readFileVar_real(file_p, 1, &time->ark_rtol);
    readFileVar_real(file_p, 1, &time->ark_atol);
    readFileVar_int( file_p, 1, &time->ark_ewtScaling);
    fscanf(file_p, "%*s");  // /.

    fscanf(file_p, "%*s");  // &species.
    readFileVar_int( file_p, 1, &pop->numSpecies);
    readFileVar_int( file_p, 1, &pop->mpiProcs);
    pop->spec = (struct species*) calloc(pop->numSpecies, sizeof(struct species));
    int *specNumMoms = (int*) calloc(pop->numSpecies, sizeof(int));
    int *specOnes    = (int*) calloc(pop->numSpecies, sizeof(int));
    int *specnDim    = (int*) calloc(pop->numSpecies, sizeof(int));
    readFileVar_int( file_p, pop->numSpecies, specNumMoms);
    int filePos = ftell(file_p);
    for (int s=0; s<pop->numSpecies; s++) {
      pop->spec[s].numMoments = specNumMoms[s];
      specOnes[s] = 1;
      specnDim[s] = nDim;
    }
    real* real_p; real** real_pp; int* int_p;
    for (int s=0; s<pop->numSpecies; s++) {
      fseek(file_p, filePos, SEEK_SET);
      real_p  =    &pop->spec[s].qCharge; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =     &pop->spec[s].muMass; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =        &pop->spec[s].tau; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =       &pop->spec[s].omSt; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =        &pop->spec[s].omd; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =      &pop->spec[s].delta; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =  &pop->spec[s].deltaPerp; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =        &pop->spec[s].eta; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_pp =      &pop->spec[s].alpha; readFileSpeciesPar_real(real_pp, file_p, s, specNumMoms);
      real_pp =         &pop->spec[s].nu; readFileSpeciesPar_real(real_pp, file_p, s, specNumMoms);
      real_p  =     &pop->spec[s].delta0; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  = &pop->spec[s].hDiffOrder; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_pp =      &pop->spec[s].hDiff; readFileSpeciesPar_real(real_pp, file_p, s, specnDim);
      real_pp =   &pop->spec[s].kDiffMin; readFileSpeciesPar_real(real_pp, file_p, s, specnDim);
      int_p   =       &pop->spec[s].icOp; readFileSpeciesPar_int(  &int_p, file_p, s, specOnes);
      real_pp =    &pop->spec[s].initAux; readFileSpeciesPar_real(real_pp, file_p, s, specnDim);
      real_p  =      &pop->spec[s].initA; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
      real_p  =     &pop->spec[s].noiseA; readFileSpeciesPar_real(&real_p, file_p, s, specOnes);
    }
    free(specnDim); free(specOnes); free(specNumMoms);
    fscanf(file_p, "%*s");  // /.
  
    fscanf(file_p, "%*s");  // &field.
    readFileVar_real(file_p, 1, &field->lambdaD);
    readFileVar_int( file_p, 1, &field->pade);
    readFileVar_int( file_p, 1, &field->icOp);
    fscanf(file_p, "%*s");  // /.

    fclose(file_p);
  }

  // Broadcast to other processes.
  MPI_Bcast(&grid->fG.Nkx     , nDim, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->fG.kxMin   , nDim, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->fG.kxMaxDyn, nDim, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&grid->mpiProcs   , nDim, mugyMPI_INT , ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&time->dt              , 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->endTime         , 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->nFrames         , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_kySplit     , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableExp, 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_fastTableImp, 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_slowTable   , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_dtFast      , 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_rtol        , 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_atol        , 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&time->ark_ewtScaling  , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);

  MPI_Bcast(&pop->numSpecies, 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&pop->mpiProcs  , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  if (myRank != ioRank) pop->spec = (struct species*) calloc(pop->numSpecies, sizeof(struct species));
  for (int s=0; s<pop->numSpecies; s++) {
    MPI_Bcast(&pop->spec[s].numMoments,                       1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
    if (myRank != ioRank) {
      pop->spec[s].alpha    = alloc_realArray(pop->spec[s].numMoments);
      pop->spec[s].nu       = alloc_realArray(pop->spec[s].numMoments);
      pop->spec[s].hDiff    = alloc_realArray(nDim);
      pop->spec[s].kDiffMin = alloc_realArray(nDim);
      pop->spec[s].initAux  = alloc_realArray(nDim);
    }
    MPI_Bcast(&pop->spec[s].muMass    ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].tau       ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].omSt      ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].omd       ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].delta     ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].deltaPerp ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].eta       ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].alpha      , pop->spec[s].numMoments, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].nu         , pop->spec[s].numMoments, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].delta0    ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].hDiffOrder,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].hDiff      ,                    nDim, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].kDiffMin   ,                    nDim, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].icOp      ,                       1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
    MPI_Bcast(pop->spec[s].initAux    ,                    nDim, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].initA     ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
    MPI_Bcast(&pop->spec[s].noiseA    ,                       1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  }

  MPI_Bcast(&field->lambdaD, 1, mugyMPI_REAL, ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->pade   , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);
  MPI_Bcast(&field->icOp   , 1, mugyMPI_INT , ioRank, MPI_COMM_WORLD);

}

void read_inputs(int argc, char *argv[], struct ioSetup *ioSet, struct grid *grid, struct timeSetup *time,
                 struct population *pop, struct fieldParameters *field) {
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

  read_inputFile(ioSet->inputFile, grid, time, pop, field);

}

void init_global_grids(struct grid *globalGrid) {
  // Set number of cells in de-aliased, aliased and real space global grids.

  /* Given user-input number of distinct dealised wavenumbers, Nkx,
     the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
     to be a power of 2, so we may need to adjust Nkx. */
  arrPrint_int(globalGrid->fG.Nkx, nDim, " User requested  NkxG=("," ) distinct wavenumbers (absolute magnitude)\n");
  globalGrid->fGa.dual.Nx[0]  = closest_power_of_two(3*(globalGrid->fG.Nkx[0]-1));
  globalGrid->fGa.dual.Nx[1]  = closest_power_of_two(3*(globalGrid->fG.Nkx[1]-1));
  globalGrid->fGa.dual.Nx[2]  = 1;
  globalGrid->fGa.dual.NxTot  = prod_int(globalGrid->fGa.dual.Nx,nDim);
  globalGrid->fGa.dual.NxyTot = prod_int(globalGrid->fGa.dual.Nx,2);

  // Number of distinct aliased (absolute) wavenumbers.
  globalGrid->fGa.Nkx[0] = globalGrid->fGa.dual.Nx[0]/2+1;
  globalGrid->fGa.Nkx[1] = globalGrid->fGa.dual.Nx[1]/2+1;
  globalGrid->fGa.Nkx[2] = globalGrid->fGa.dual.Nx[2]/2+1;
  // Length of aliased arrays along kx and ky.
  globalGrid->fGa.Nekx[0]  = globalGrid->fGa.dual.Nx[0];
  globalGrid->fGa.Nekx[1]  = globalGrid->fGa.Nkx[1];
  globalGrid->fGa.Nekx[2]  = globalGrid->fGa.Nkx[2];
  globalGrid->fGa.NekxTot  = prod_int(globalGrid->fGa.Nekx,nDim);
  globalGrid->fGa.NekxyTot = prod_int(globalGrid->fGa.Nekx,2);

  // Recompute the number of distinct de-aliased (absolute) wavenumbers.
  globalGrid->fG.Nkx[0] = 2*(globalGrid->fGa.Nkx[0]-1)/3+1;
  globalGrid->fG.Nkx[1] = 2*(globalGrid->fGa.Nkx[1]-1)/3+1;
  globalGrid->fG.Nkx[2] = 1;
  // Length of de-aliased arrays along kx and ky.
  globalGrid->fG.Nekx[0]  = 2*(globalGrid->fG.Nkx[0]-1)+1;
  globalGrid->fG.Nekx[1]  = globalGrid->fG.Nkx[1];
  globalGrid->fG.Nekx[2]  = globalGrid->fG.Nkx[2];
  globalGrid->fG.NekxTot  = prod_int(globalGrid->fG.Nekx,nDim);
  globalGrid->fG.NekxyTot = prod_int(globalGrid->fG.Nekx,2);

  // Number of cells in de-aliased real-space.
  globalGrid->fG.dual.Nx[0]  = 2*(globalGrid->fG.Nkx[0]-1)+1;
  globalGrid->fG.dual.Nx[1]  = 2*(globalGrid->fG.Nkx[1]-1);
  globalGrid->fG.dual.Nx[2]  = 1;
  globalGrid->fG.dual.NxTot  = prod_int(globalGrid->fG.dual.Nx,nDim);
  globalGrid->fG.dual.NxyTot = prod_int(globalGrid->fG.dual.Nx,2);

  real Lx[nDim] = {2.0*M_PI/globalGrid->fG.kxMin[0], 2.0*M_PI/globalGrid->fG.kxMin[1], 2.0*M_PI/globalGrid->fG.kxMin[2]};

  // Length of dealised and aliased real-space cell.
  for (int d=0; d<nDim; d++) {
    globalGrid->fG.dual.dx[d]  = Lx[d]/fmax(1.,(real)(globalGrid->fG.dual.Nx[d]-globalGrid->fG.dual.Nx[d] % 2));
    globalGrid->fGa.dual.dx[d] = Lx[d]/fmax(1.,(real)(globalGrid->fGa.dual.Nx[d]-globalGrid->fGa.dual.Nx[d] % 2));
  }

  // Global de-aliased real-space grids
  globalGrid->fG.dual.x = alloc_realArray(sum_int(globalGrid->fG.dual.Nx, nDim));
  real *dx = globalGrid->fG.dual.dx;
  int *Nx = globalGrid->fG.dual.Nx; 
  int xOff = 0;
  for (int d=0; d<nDim; d++) {
    for (int i=0; i<Nx[d]; i++)
      globalGrid->fG.dual.x[i+xOff] = (real)(i)*dx[d]+(real)(1-Nx[d] % 2-1)*0.5*dx[d]-0.5*Lx[d];
    globalGrid->fG.dual.xMin[d] = globalGrid->fG.dual.x[0+xOff];
    globalGrid->fG.dual.xMax[d] = globalGrid->fG.dual.x[Nx[d]-1+xOff];
    xOff += Nx[d];
  }
  // Global aliased real-space grids (may not be needed).
  globalGrid->fGa.dual.x = alloc_realArray(sum_int(globalGrid->fGa.dual.Nx, 3));
  real *dxa = globalGrid->fGa.dual.dx;
  int *Nxa = globalGrid->fGa.dual.Nx; 
  int xaOff = 0;
  for (int d=0; d<nDim; d++) {
    for (int i=0; i<Nxa[d]; i++)
      globalGrid->fGa.dual.x[i+xaOff] = (real)(i)*dxa[d]+(real)(1-Nxa[d] % 2-1)*0.5*dxa[d]-0.5*Lx[d];
    globalGrid->fGa.dual.xMin[d] = globalGrid->fGa.dual.x[0+xaOff];
    globalGrid->fGa.dual.xMax[d] = globalGrid->fGa.dual.x[Nxa[d]-1+xaOff];
    xaOff += Nxa[d];
  }

  // Global dealiased k-space grids.
  for (int d=0; d<nDim; d++) globalGrid->fGa.kxMin[d] = globalGrid->fG.kxMin[d];
  globalGrid->fG.kx  = alloc_realArray(sum_int(globalGrid->fG.Nekx, 3));
  real *kxMin = globalGrid->fG.kxMin;
  int *Nkx = globalGrid->fG.Nkx; 
  int kxOff = 0;
  for (int d=0; d<nDim; d++) {
    for (int i=0; i<Nkx[d]; i++)
      globalGrid->fG.kx[i+kxOff] = (real)(i)*kxMin[d];
    kxOff += globalGrid->fG.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (int i=Nkx[0]; i<globalGrid->fG.Nekx[0]; i++)
    globalGrid->fG.kx[i] = -(real)(Nkx[0]-1-(i-Nkx[0]))*kxMin[0];

  // Global aliased k-space grids.
  globalGrid->fGa.kx = alloc_realArray(sum_int(globalGrid->fGa.Nekx, 3));
  real *kxaMin = globalGrid->fGa.kxMin;
  int *Nkxa = globalGrid->fGa.Nkx; 
  int kxaOff = 0;
  for (int d=0; d<nDim; d++) {
    for (int i=0; i<Nkxa[d]; i++)
      globalGrid->fGa.kx[i+kxaOff] = (real)(i)*kxaMin[d];
    kxaOff += globalGrid->fGa.Nekx[d];
  }
  // Negative kx modes in increasing order.
  for (int i=Nkxa[0]; i<globalGrid->fGa.Nekx[0]; i++)
    globalGrid->fGa.kx[i] = -(real)(Nkxa[0]-1-(i-Nkxa[0]))*kxaMin[0];

  r0printf("\n Proceeding with :\n");
  arrPrint_int(globalGrid->fG.Nkx,      nDim, " Number of distinct de-aliased absolute wavenumbers: NkxG   =", "\n");
  arrPrint_int(globalGrid->fG.Nekx,     nDim, " Length of de-aliased k-space arrays:                NekxG  =", "\n");
  arrPrint_int(globalGrid->fGa.Nkx,     nDim, " Number of distinct aliased absolute wavenumbers:    NkxaG  =", "\n");
  arrPrint_int(globalGrid->fGa.Nekx,    nDim, " Length of aliased k-space arrays:                   NekxaG =", "\n");
  arrPrint_int(globalGrid->fGa.dual.Nx, nDim, " Number of aliased real space cells:                 NxaG   =", "\n");
  arrPrint_int(globalGrid->fG.dual.Nx,  nDim, " Number of de-aliased real space cells:              NxG    =", "\n");

  arrPrint_real(globalGrid->fG.kxMin,    nDim, " Minimum absolute magnitude of wavenumbers: kxMin    =", "\n");
  arrPrint_real(globalGrid->fG.kxMaxDyn, nDim, " Largest wavenumbers evolved:               kxMaxDyn =", "\n");

}

void allocate_fields(struct grid localGrid, struct population localPop) {
  // Allocate various fields needed.
  resource onResource = hostOnly;
  alloc_fourierMoments( localGrid.fG, localPop, onResource, &momk);
}

void set_initialCondition(struct grid localGrid, struct population localPop) {
  // Impose the initial conditions on the moments and thoe potential.

  int sOff = 0;
  real initA = 2.5e-2;
  real initAux = -3.;
  real *kxMin = &localGrid.fG.kxMin[0];
  for (int s=0; s<localPop.numSpecies; s++) {
    fourier *den_p  = getMoment_fourier(localGrid.fG, localPop, s, denIdx, momk.ho);  // Get density of species s.
    fourier *temp_p = getMoment_fourier(localGrid.fG, localPop, s, tempIdx, momk.ho);  // Get temperature of species s.
    for (int linIdx=0; linIdx<localGrid.fG.NekxTot; linIdx++) {
      int kxIdx[nDim];
      lin2sub_fourier(&kxIdx[0], linIdx, localGrid.fG);  // Convert linear index to multidimensional kx index.
      real kx[nDim];
      get_kx(&kx[0], kxIdx, localGrid.fG);

      // Set density to a power-law in k-space.
      den_p++;
      den_p[0] = initA*(pow((kxMin[0]+abs(kx[0]))/kxMin[0],initAux))
                      *(pow((kxMin[1]+abs(kx[1]))/kxMin[1],initAux));

      // Set the initial temperature (fluctuations) to zero.
      temp_p++;
      temp_p[0] = 0.;
    };

  }
}

void free_fields() {
  // Deallocate fields.
  resource onResource = hostOnly;
  free_fourierMoments(&momk, onResource);
}

void free_grid(struct grid *grid) {
  // Deallocate memory used by grids.
  free(grid->fG.dual.x);
  free(grid->fG.kx);
  free(grid->fGa.dual.x);
  free(grid->fGa.kx);
}

void free_population(struct population *pop) {
  // Deallocate memory used in species struct.
  for (int s=0; s<pop->numSpecies; s++) {
    free(pop->spec[s].alpha);
    free(pop->spec[s].nu);
    free(pop->spec[s].hDiff);
    free(pop->spec[s].kDiffMin);
    free(pop->spec[s].initAux);
  }
  free(pop->spec);
}
