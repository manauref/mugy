#!/bin/bash


# determine machine "name" based on 'HOST' environment variable, or
# if that is unset, then via the shell 'hostname' command
HOST=${HOST:-`hostname`}

# set hostname-specific defaults
if [[ $HOST == "mfrancis-lt" ]]
then

  #[ Mana's local.
  export CCCOMPILER=mpicc
  export CCFLAGS=
  export MPI_DIR=$HOME/Documents/multiscale/code/openmpi4.1.1-intel19.1//
  export CCCOMPILER=$MPI_DIR/bin/mpicc
  export MPI_INC=$MPI_DIR/include
  export MPI_LIB=$MPI_DIR/lib
  export CUDA_DIR=
  export CUDA_INC=
  export CUDA_LIB=
  export ADIOS_DIR=$HOME/Documents/multiscale/code/adios2-openmpi4.1.1-intel19.1/
  export ADIOS_INC=`$ADIOS_DIR/bin/adios2-config --c-flags`
  export ADIOS_LIB=`$ADIOS_DIR/bin/adios2-config --c-libs`
  export SUNDIALS_DIR=$HOME/multiscale/code/BUILD.Release/INSTALL.Release/
  export SUNDIALS_MOD=$SUNDIALS_DIR/fortran
  export SUNDIALS_LIB=$SUNDIALS_DIR/lib

elif [[ "$HOST" == *"stellar"* ]]
then

  #[ Princeton's Stellar.
  module load cudatoolkit/11.1
  module load openmpi/cuda-11.1/gcc/4.1.1
  module load anaconda3/2020.11
  module load fftw/nvhpc-21.5/openmpi-4.1.1/3.3.9
  export MPI_DIR=$MPI_HOME
  export CCCOMPILER=$MPI_DIR/bin/mpicc
  export NVCCCOMPILER=$CUDA_HOME/bin/nvcc
  export NVCCFLAGS=-arch=sm_80
  export CCFLAGS=-lm
  export MPI_INC=$MPI_DIR/include
  export MPI_LIB=$MPI_DIR/lib64
  export FFTW_INC=$FFTW3DIR/include
  export FFTW_LIB=$FFTW3DIR/lib64
  export CUDA_DIR=$CUDA_HOME
  export CUDA_INC=$CUDA_DIR/include
  export CUDA_LIB=$CUDA_DIR/lib64
  export ADIOS_DIR=$HOME/multiscale/code/adios2-openmpi-cuda-11.1-gcc-4.1.1/
  export ADIOS_INC=`$ADIOS_DIR/bin/adios2-config --c-flags`
  export ADIOS_LIB=`$ADIOS_DIR/bin/adios2-config --c-libs`
  export SUNDIALS_DIR=$HOME/multiscale/code/BUILD.Release/INSTALL.Release/
  export SUNDIALS_MOD=$SUNDIALS_DIR/fortran
  export SUNDIALS_LIB=$SUNDIALS_DIR/lib

else

  echo ''
#  #[ Darin Ernst's computer:
#  export FCOMPILER=gfortran-mp-6
#  export FFTW_DIR=/opt/local

#  #[ Mana's TheFarm:
#  export FCOMPILER=/usr/local/Cellar/gcc/8.2.0/bin/gfortran-8
#  export FFTW_DIR=/usr/local/Cellar/fftw/3.3.8/

#  #[ Plankton:
#  export FCOMPILER=/opt/intel/compilers_and_libraries_2017.7.259/linux/bin/intel64/ifort
#  export FFTW_DIR=/usr/local/

#  #[ NERSC's Cori:
#  module load intel
#  module load fftw/3.3.8
#  module load adios/1.13.1
#  export FCOMPILER=mpiifort
#  export ADIOS_FREADLIB=$ADIOSREAD_FLIB

#  #[ Steve Leak's setup for NERSC's Cori (temporary):
#  module load adios cray-fftw
#  export ADIOS_FREADLIB="-L/global/common/sw/cray/cnl7/haswell/adios/1.13.1/intel/19.0.3.199/wuurn4w/lib -L/global/common/sw/cray/cnl7/haswell/zfp/0.5.0/intel/19.0.3.199/cle2mg7/lib -L/global/common/sw/cray/cnl7/haswell/c-blosc/1.16.3/intel/19.0.3.199/tn4uoyl/lib -L/global/common/sw/cray/cnl7/haswell/zstd/1.4.0/intel/19.0.3.199/k2svlvb/lib -L/global/common/sw/cray/cnl7/haswell/snappy/1.1.7/intel/19.0.3.199/vlbmtjf/lib -L/global/common/sw/cray/cnl7/haswell/lz4/1.9.0/intel/19.0.3.199/qf6mmqo/lib -L/global/common/sw/cray/cnl7/haswell/sz/1.4.12.3/intel/19.0.3.199/pubgujg/lib -L/global/common/sw/cray/cnl7/haswell/snappy/1.1.7/intel/19.1.2.254/6e4sd6n/lib -L/global/common/sw/cray/cnl7/haswell/bzip2/1.0.6/intel/19.0.3.199/p7h55tt/lib/ -ladiosreadf -ladiosf -lblosc -lsnappy -lzfp -lbz2 -llz4 -lbz2 -lSZ -lzstd -cxxlib"
#  export ADIOS_FLIB="-L/global/common/sw/cray/cnl7/haswell/adios/1.13.1/intel/19.0.3.199/wuurn4w/lib -ladiosf"
#  export FCOMPILER=ftn

fi
