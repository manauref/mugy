/* mugy: mpi_tools

   MPI-related operations.
*/

#include "mpi_tools.h"

int myRank, totNumProcs;  // Rank of this process & total number of processes.
int numProcs[nDim+1];     // Processes along x,y,z and species.
MPI_Comm cartComm;        // Cartesian communicator.
int cartRank;             // Rank in cartCOMM.
MPI_Comm *sub1dComm;      // 1D subcommunicators along each direction.
int sub1dRank[nDim+1];    // ID (rank) in the 1D xpec,Z,X,Y subcommunicators.
MPI_Comm *sComm, *zComm, *xComm, *yComm;  // Pointers for 1D comms.
int sRank, zRank, xRank, yRank;  // Pointers for 1D rank IDs.

void init_mpi(int argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  // Rank of this MPI process.
  MPI_Comm_size(MPI_COMM_WORLD, &totNumProcs);  // Number of MPI processes.
}

void init_comms(struct grid grid, struct speciesParameters spec) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_int(grid.mpiProcs,nDim)*spec.mpiProcs != totNumProcs) {
    printf(" Number of MPI processes in input file (%d %d) differs from that in mpirun (%d).\n",
           prod_int(grid.mpiProcs,nDim),spec.mpiProcs, totNumProcs);
    abortSimulation(" Terminating...\n");
  }

  // Number of MPI processes along spec,Z,X,Y.
  numProcs[0] = spec.mpiProcs;
  numProcs[1] = grid.mpiProcs[2];
  numProcs[2] = grid.mpiProcs[0];
  numProcs[3] = grid.mpiProcs[1];
  arrPrint_int(numProcs,nDim+1, " MPI processes along spec,Z,X,Y: ", "\n");
  r0printf("\n");

  // MPI_CART boundary conditions. True=periodic.
  int cartCommBCs[nDim+1] = {true, true, true, false};
  // Let MPI assign arbitrary ranks.
  int reorder = true;

  // Create a 4D Cartesian communicator (3D space + species).
  MPI_Cart_create(MPI_COMM_WORLD, nDim+1, numProcs, cartCommBCs, reorder, &cartComm);

  // Find my coordinate parameters in the Cartesial topology.
  MPI_Comm_rank(cartComm, &cartRank);
  int cartCoords[nDim+1];
  MPI_Cart_coords(cartComm, myRank, nDim+1, &cartCoords[0]);

  // spec,Z,X,Y SUBCOMMUNICATORS.
  sub1dComm = (MPI_Comm *) calloc(nDim+1, sizeof(MPI_Comm));
  for (int d=0; d<nDim+1; d++) {
    int remain[nDim+1] = {false,false,false,false};
    remain[d] = true;     // Keep lc direction.
    MPI_Cart_sub(cartComm, remain, &sub1dComm[d]);
    MPI_Comm_rank(sub1dComm[d], &sub1dRank[d]);
  }
  sComm = &sub1dComm[0];  zComm = &sub1dComm[1];
  xComm = &sub1dComm[2];  yComm = &sub1dComm[3];
  sRank = sub1dRank[0];  zRank = sub1dRank[1];
  xRank = sub1dRank[2];  yRank = sub1dRank[3];

}

void distribute1dDOFs(const int procs, const int procID, const int globalDOFs, int *localDOFs, int *firstDOF) {
  /* Distribute 1D degrees of freedom amongst processes. Give every process at
     least the floor of the division of the degrees of freedom by the number
     of processors. Then distribute the remaining elements one at a time.
       procs:  number of processes
       procID: this process' rank.
       globalDOFs: number of degrees of freedom in global domain.
       localDOFs: number of degrees of freedom in this process.
       firstDOF: first degree of freedom from the global DOFs in this process.
   */
  *localDOFs = floor((real)(globalDOFs)/(real)(procs));
  int remainingDOFs = globalDOFs - procs*(*localDOFs);
  *firstDOF = procID*(*localDOFs);

  while (remainingDOFs > 0) {
    for (int p=0; p<procs; p++) {
      remainingDOFs += -1;
      if (p == procID) {
        *localDOFs += 1;
      } else if (p < procID) {
        *firstDOF += 1;
      }
      if (remainingDOFs == 0) break;
    }
  };

}

void allocAndCopyVar_int(int **var, int *src, const int numE) {
  // Allocate int array var and copy numE elements to it from src.
  *var = alloc_intArray(numE);
  memcpy(&(*var)[0], src, numE*sizeof(int));
}

void allocAndCopyVar_real(real **var, real *src, const int numE) {
  // Allocate real array var and copy numE elements to it from src.
  *var = alloc_realArray(numE);
  memcpy(&(*var)[0], src, numE*sizeof(real));
}

void distributeDOFs(struct grid globalGrid, struct speciesParameters globalSpec, struct grid *localGrid, struct speciesParameters *localSpec) {
  // Distribute s,Z,X,Y amongst MPI processes.
  
  // Distribute the species.
  distribute1dDOFs(globalSpec.mpiProcs, sRank, globalSpec.numSpecies, &localSpec->numSpecies, &localSpec->globalOff);
  localSpec->numMoments = globalSpec.numMoments;
  allocAndCopyVar_real(&localSpec->qCharge   ,globalSpec.qCharge   +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->muMass    ,globalSpec.muMass    +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->tau       ,globalSpec.tau       +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->omSt      ,globalSpec.omSt      +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->omd       ,globalSpec.omd       +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->delta     ,globalSpec.delta     +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->deltaPerp ,globalSpec.deltaPerp +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->eta       ,globalSpec.eta       +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->alpha     ,globalSpec.alpha     +localSpec->numMoments*localSpec->globalOff,localSpec->numMoments*localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->nu        ,globalSpec.nu        +localSpec->numMoments*localSpec->globalOff,localSpec->numMoments*localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->delta0    ,globalSpec.delta0    +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->hDiffOrder,globalSpec.hDiffOrder+                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->hDiff     ,globalSpec.hDiff     +                 nDim*localSpec->globalOff,                 nDim*localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->kDiffMin  ,globalSpec.kDiffMin  +                 nDim*localSpec->globalOff,                 nDim*localSpec->numSpecies);
  allocAndCopyVar_int(&localSpec->icOp      ,globalSpec.icOp      +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->initAux   ,globalSpec.initAux   +                 nDim*localSpec->globalOff,                 nDim*localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->initA     ,globalSpec.initA     +                      localSpec->globalOff,                      localSpec->numSpecies);
  allocAndCopyVar_real(&localSpec->noiseA    ,globalSpec.noiseA    +                      localSpec->globalOff,                      localSpec->numSpecies);

  // Distribute the real-space and Fourier-space points. 
  int dualD[nDim] = {2,3,1};
  for (int d=0; d<nDim; d++) {
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[dualD[d]], globalGrid.fG.Nekx[d],
                     &localGrid->fG.Nekx[d], &localGrid->fG.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[dualD[d]], globalGrid.fG.dual.Nx[d],
                     &localGrid->fG.dual.Nx[d], &localGrid->fG.dual.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[dualD[d]], globalGrid.fGa.Nekx[d],
                     &localGrid->fGa.Nekx[d], &localGrid->fGa.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[dualD[d]], globalGrid.fGa.dual.Nx[d],
                     &localGrid->fGa.dual.Nx[d], &localGrid->fGa.dual.globalOff[d]);
  }

  // Create local real-space and Fourier-space coordinate arrays.
  localGrid->fG.kx     = alloc_realArray(sum_int(localGrid->fG.Nekx, nDim));
  localGrid->fG.dual.x = alloc_realArray(sum_int(localGrid->fG.dual.Nx, nDim));
  localGrid->fGa.kx     = alloc_realArray(sum_int(localGrid->fGa.Nekx, nDim));
  localGrid->fGa.dual.x = alloc_realArray(sum_int(localGrid->fGa.dual.Nx, nDim));
  for (int d=0; d<nDim; d++) {
    memcpy(getArray_real(localGrid->fG.kx,localGrid->fG.Nekx,d),
           getArray_real(globalGrid.fG.kx,globalGrid.fG.Nekx,d)+localGrid->fG.globalOff[d], localGrid->fG.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fG.dual.x,localGrid->fG.dual.Nx,d),
           getArray_real(globalGrid.fG.dual.x,globalGrid.fG.dual.Nx,d)+localGrid->fG.dual.globalOff[d], localGrid->fG.dual.Nx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.kx,localGrid->fGa.Nekx,d),
           getArray_real(globalGrid.fGa.kx,globalGrid.fGa.Nekx,d)+localGrid->fGa.globalOff[d], localGrid->fGa.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.dual.x,localGrid->fGa.dual.Nx,d),
           getArray_real(globalGrid.fGa.dual.x,globalGrid.fGa.dual.Nx,d)+localGrid->fGa.dual.globalOff[d], localGrid->fGa.dual.Nx[d]*sizeof(real));
  }

//  int dO = 2;
//  arrPrint_real(getArray_real(globalGrid.fG.kx,globalGrid.fG.Nekx,dO),globalGrid.fG.Nekx[dO], " kxG = ", "\n");
//  r0printf("\n");
//  int maxCharacters = 9999;
//  char strOut[maxCharacters];
//  arr2str_real(strOut, getArray_real(localGrid->fG.kx,localGrid->fG.Nekx,dO), localGrid->fG.Nekx[dO], "kxL = ", "\n");
//  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.
//  if (myRank==0) { printf("r0 : \n"); printf(strOut); }
//  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.
//  if (myRank==1) { printf("r1 : \n"); printf(strOut); }
//  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

}

void terminate_mpi() {
  // Close MPI up, thoroughly.
  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  for (int d=0; d<nDim+1; d++) MPI_Comm_free(&sub1dComm[d]);
  free(sub1dComm);

  MPI_Comm_free(&cartComm);

  MPI_Finalize();
}
