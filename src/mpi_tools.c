/* mugy: mpi_tools

   MPI-related operations.
*/

#include "mpi_tools.h"

mint myRank, totNumProcs;  // Rank of this process & total number of processes.
mint numProcs[nDim+1];     // Processes along x,y,z and species.
MPI_Comm cartComm;        // Cartesian communicator.
mint cartRank;             // Rank in cartCOMM.
MPI_Comm *sub1dComm;      // 1D subcommunicators along each direction.
mint sub1dRank[nDim+1];    // ID (rank) in the 1D xpec,Z,X,Y subcommunicators.
MPI_Comm *sComm, *zComm, *xComm, *yComm;  // Pointers for 1D comms.
mint sRank, zRank, xRank, yRank;  // Pointers for 1D rank IDs.

void init_mpi(mint argc, char *argv[]) {
  // Initialize MPI, get rank of this process and total number of processes.
  MPI_Init(&argc, &argv);  // Initialize MPI.
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  // Rank of this MPI process.
  MPI_Comm_size(MPI_COMM_WORLD, &totNumProcs);  // Number of MPI processes.
}

void init_comms(struct grid grid, struct population pop) {
  // Initialize the various sub-communicators needed.
  
  // Check the number of MPI processes is correct.
  if (prod_mint(grid.mpiProcs,nDim)*pop.mpiProcs != totNumProcs) {
    printf(" Number of MPI processes in input file (%d) differs from that in mpirun (%d).\n",
           prod_mint(grid.mpiProcs,nDim)*pop.mpiProcs, totNumProcs);
    abortSimulation(" Terminating...\n");
  }

  // Number of MPI processes along X,Y,Z,s.
  numProcs[0] = grid.mpiProcs[0];
  numProcs[1] = grid.mpiProcs[1];
  numProcs[2] = grid.mpiProcs[2];
  numProcs[3] = pop.mpiProcs;
  arrPrint_mint(numProcs,nDim+1, " MPI processes along X,Y,Z,s: ", "\n");
  r0printf("\n");

  // MPI_CART boundary conditions. True=periodic.
  mint cartCommBCs[nDim+1] = {true, true, true, false};
  // Let MPI assign arbitrary ranks.
  mint reorder = true;

  // Create a 4D Cartesian communicator (3D space + species).
  // Re-organize the dimensions in s,Z,X,Y order.
  mint commOrg[nDim+1] = {2,3,1,0};
  mint numProcsOrg[nDim+1], cartCommBCsOrg[nDim+1];
  for (mint d=0; d<nDim+1; d++) {
    numProcsOrg[commOrg[d]] = numProcs[d];
    cartCommBCsOrg[commOrg[d]] = cartCommBCs[d];
  }
  MPI_Cart_create(MPI_COMM_WORLD, nDim+1, numProcsOrg, cartCommBCsOrg, reorder, &cartComm);

  // Find my coordinate parameters in the Cartesial topology.
  MPI_Comm_rank(cartComm, &cartRank);
  mint cartCoords[nDim+1];
  MPI_Cart_coords(cartComm, myRank, nDim+1, &cartCoords[0]);

  // 1D subcommunicators. Organize them in s,Z,X,Y, order so processes of a
  // single X-Y plane of a single species are contiguous.
  sub1dComm = (MPI_Comm *) calloc(nDim+1, sizeof(MPI_Comm));
  for (mint d=0; d<nDim+1; d++) {
    mint remain[nDim+1] = {false,false,false,false};
    remain[commOrg[d]] = true;
    MPI_Cart_sub(cartComm, remain, &sub1dComm[d]);
    MPI_Comm_rank(sub1dComm[d], &sub1dRank[d]);
  }
  xComm = &sub1dComm[0];  yComm = &sub1dComm[1];
  zComm = &sub1dComm[2];  sComm = &sub1dComm[3];
  xRank = sub1dRank[0];  yRank = sub1dRank[1];
  zRank = sub1dRank[2];  sRank = sub1dRank[3];

}

void distribute1dDOFs(const mint procs, const mint procID, const mint globalDOFs, mint *localDOFs, mint *firstDOF) {
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
  mint remainingDOFs = globalDOFs - procs*(*localDOFs);
  *firstDOF = procID*(*localDOFs);

  while (remainingDOFs > 0) {
    for (mint p=0; p<procs; p++) {
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

void distributeDOFs(struct grid globalGrid, struct population globalPop, struct grid *localGrid, struct population *localPop) {
  // Distribute s,Z,X,Y amongst MPI processes.
  
  // Distribute the species.
  distribute1dDOFs(globalPop.mpiProcs, sRank, globalPop.numSpecies, &localPop->numSpecies, &localPop->globalSpecOff);
  localPop->globalMomOff = 0;
  for (mint s=0; s<localPop->globalSpecOff; s++) localPop->globalMomOff += globalPop.spec[s].numMoments;
  localPop->spec = (struct species*) calloc(localPop->numSpecies, sizeof(struct species));
  for (mint s=0; s<localPop->numSpecies; s++) {
    localPop->spec[s].numMoments = globalPop.spec[s+localPop->globalSpecOff].numMoments;

    localPop->spec[s].alpha      = alloc_realArray(localPop->spec[s].numMoments);
    localPop->spec[s].nu         = alloc_realArray(localPop->spec[s].numMoments);
    localPop->spec[s].hDiffOrder = alloc_realArray(nDim);
    localPop->spec[s].hDiff      = alloc_realArray(nDim);
    localPop->spec[s].kDiffMin   = alloc_realArray(nDim);
    localPop->spec[s].initAux    = alloc_realArray(nDim);

    localPop->spec[s].qCharge    = globalPop.spec[s+localPop->globalSpecOff].qCharge   ;
    localPop->spec[s].muMass     = globalPop.spec[s+localPop->globalSpecOff].muMass    ;
    localPop->spec[s].tau        = globalPop.spec[s+localPop->globalSpecOff].tau       ;
    localPop->spec[s].omSt       = globalPop.spec[s+localPop->globalSpecOff].omSt      ;
    localPop->spec[s].omd        = globalPop.spec[s+localPop->globalSpecOff].omd       ;
    localPop->spec[s].delta      = globalPop.spec[s+localPop->globalSpecOff].delta     ;
    localPop->spec[s].deltaPerp  = globalPop.spec[s+localPop->globalSpecOff].deltaPerp ;
    localPop->spec[s].eta        = globalPop.spec[s+localPop->globalSpecOff].eta       ;
    memcpy(localPop->spec[s].alpha, globalPop.spec[s+localPop->globalSpecOff].alpha, localPop->spec[s].numMoments*sizeof(real));
    memcpy(localPop->spec[s].nu   , globalPop.spec[s+localPop->globalSpecOff].nu   , localPop->spec[s].numMoments*sizeof(real));
    localPop->spec[s].delta0     = globalPop.spec[s+localPop->globalSpecOff].delta0    ;
    memcpy(localPop->spec[s].hDiffOrder, globalPop.spec[s+localPop->globalSpecOff].hDiffOrder, nDim*sizeof(real));
    memcpy(localPop->spec[s].hDiff     , globalPop.spec[s+localPop->globalSpecOff].hDiff     , nDim*sizeof(real));
    memcpy(localPop->spec[s].kDiffMin  , globalPop.spec[s+localPop->globalSpecOff].kDiffMin  , nDim*sizeof(real));
    localPop->spec[s].icOp       = globalPop.spec[s+localPop->globalSpecOff].icOp      ;
    memcpy(localPop->spec[s].initAux, globalPop.spec[s+localPop->globalSpecOff].initAux, nDim*sizeof(real));
    localPop->spec[s].initA      = globalPop.spec[s+localPop->globalSpecOff].initA     ;
    localPop->spec[s].noiseA     = globalPop.spec[s+localPop->globalSpecOff].noiseA    ;
  }
  // Set the total number of moments.
  localPop->numMomentsTot = 0;
  for (int s=0; s<localPop->numSpecies; s++) localPop->numMomentsTot += localPop->spec[s].numMoments;

  // Distribute the real-space and Fourier-space points. 
  for (mint d=0; d<nDim; d++) {
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[d], globalGrid.fG.Nekx[d],
                     &localGrid->fG.Nekx[d], &localGrid->fG.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[d], globalGrid.fG.dual.Nx[d],
                     &localGrid->fG.dual.Nx[d], &localGrid->fG.dual.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[d], globalGrid.fGa.Nekx[d],
                     &localGrid->fGa.Nekx[d], &localGrid->fGa.globalOff[d]);
    distribute1dDOFs(globalGrid.mpiProcs[d], sub1dRank[d], globalGrid.fGa.dual.Nx[d],
                     &localGrid->fGa.dual.Nx[d], &localGrid->fGa.dual.globalOff[d]);
  }
  localGrid->fG.NekxTot      = prod_mint(localGrid->fG.Nekx,nDim);
  localGrid->fG.dual.NxTot   = prod_mint(localGrid->fG.dual.Nx,nDim);
  localGrid->fG.NekxyTot     = prod_mint(localGrid->fG.Nekx,2);
  localGrid->fG.dual.NxyTot  = prod_mint(localGrid->fG.dual.Nx,2);
  localGrid->fGa.NekxTot     = prod_mint(localGrid->fGa.Nekx,nDim);
  localGrid->fGa.dual.NxTot  = prod_mint(localGrid->fGa.dual.Nx,nDim);
  localGrid->fGa.NekxyTot    = prod_mint(localGrid->fGa.Nekx,2);
  localGrid->fGa.dual.NxyTot = prod_mint(localGrid->fGa.dual.Nx,2);

  // Create local real-space and Fourier-space coordinate arrays.
  localGrid->fG.kx     = alloc_realArray(sum_mint(localGrid->fG.Nekx, nDim));
  localGrid->fG.dual.x = alloc_realArray(sum_mint(localGrid->fG.dual.Nx, nDim));
  localGrid->fGa.kx     = alloc_realArray(sum_mint(localGrid->fGa.Nekx, nDim));
  localGrid->fGa.dual.x = alloc_realArray(sum_mint(localGrid->fGa.dual.Nx, nDim));
  for (mint d=0; d<nDim; d++) {
    memcpy(getArray_real(localGrid->fG.kx,localGrid->fG.Nekx,d),
           getArray_real(globalGrid.fG.kx,globalGrid.fG.Nekx,d)+localGrid->fG.globalOff[d], localGrid->fG.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fG.dual.x,localGrid->fG.dual.Nx,d),
           getArray_real(globalGrid.fG.dual.x,globalGrid.fG.dual.Nx,d)+localGrid->fG.dual.globalOff[d], localGrid->fG.dual.Nx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.kx,localGrid->fGa.Nekx,d),
           getArray_real(globalGrid.fGa.kx,globalGrid.fGa.Nekx,d)+localGrid->fGa.globalOff[d], localGrid->fGa.Nekx[d]*sizeof(real));
    memcpy(getArray_real(localGrid->fGa.dual.x,localGrid->fGa.dual.Nx,d),
           getArray_real(globalGrid.fGa.dual.x,globalGrid.fGa.dual.Nx,d)+localGrid->fGa.dual.globalOff[d], localGrid->fGa.dual.Nx[d]*sizeof(real));
  }

  // Copy global scalars into local grids.
  for (mint d=0; d<nDim; d++) {
    localGrid->fG.kxMin[d]    = globalGrid.fG.kxMin[d];
    localGrid->fGa.kxMin[d]   = globalGrid.fGa.kxMin[d];
    localGrid->fG.dual.dx[d]  = globalGrid.fG.dual.dx[d];
    localGrid->fGa.dual.dx[d] = globalGrid.fGa.dual.dx[d];
  }

}

void terminate_mpi() {
  // Close MPI up, thoroughly.
  MPI_Barrier(MPI_COMM_WORLD); // To avoid premature deallocations.

  for (mint d=0; d<nDim+1; d++) MPI_Comm_free(&sub1dComm[d]);
  free(sub1dComm);

  MPI_Comm_free(&cartComm);

  MPI_Finalize();
}
