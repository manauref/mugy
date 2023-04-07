#[
#[ Makefile for mugy
#[
#[
#[ Users need to set the following environment variables:
#[
#[      CCOMPILER,
#[      FFTW_INC, FFTW_LIB,
#[      ADIOS_DIR, ADIOS_INC, ADIOS_FLAGS,
#[      NVCCCOMPILER, NVCCFLAGS, CUDA_INC, CUDA_LIB, 
#[
#[ This is done in setup.sh, which is invoked using ". setup.sh"
#[ Below, defaults are used if the variables are not set.
#[
#[
#[ 'make depend' uses makedepend to automatically generate dependencies
#[               (dependencies are added to end of Makefile)
#[ 'make'        build executable file 'mugy.x'
#[ 'make clean'  removes all .o and executable files
#[


#.Extract macros from header file.
USE_SINGLE_PRECISION := $(shell grep -e '\#define USE_SINGLE_PRECISION' ../include/mh_userFLAGS.h | sed 's/\#define USE_SINGLE_PRECISION //' | sed -r 's/\s+//g' )

#[ Define the C compiler to use.
CC = ${CCCOMPILER}

#[ CUDA flags and variables.
NVCC_RESULT := $(shell which nvcc 2> NULL)
NVCC_TEST := $(notdir $(NVCC_RESULT))
ifeq ($(NVCC_TEST),nvcc)
    USE_GPU = true
    CUDA_INCL = -I${CUDA_INC}
    CUDA_LIBS = -L${CUDA_LIB}
    CUDAFLAGS = -DUSE_GPU=$(USE_GPU) ${NVCCFLAGS}
    NVCC = ${NVCCCOMPILER}
    CUDART = -lcudart -lcudadevrt
    CFLAGS_DEV = -DUSE_GPU=$(USE_GPU) 
else
    USE_GPU = false
    CUDA_INCL =
    CUDA_LIBS =
    CUDAFLAGS =
    NVCC =
    CUDART =
    CFLAGS_DEV =
endif

#[ Use a different FFTW library depending on precision.
ifeq ($(USE_SINGLE_PRECISION), true)
    FFTW_LIB_FLAGS = -lfftw3f_mpi -lfftw3f -lm
else
    FFTW_LIB_FLAGS = -lfftw3_mpi -lfftw3 -lm
endif

#[ Define any compile-time flags.
#CFLAGS = -Wall -std=c11 -lm $(ADIOS_FLAGS)
#CFLAGS = -Wall -std=c11 -lm
CFLAGS = $(CFLAGS_DEV) -Wall -std=c11 ${CCFLAGS}

#[ Define any directories containing header files.
INCLUDES = -I../include -I${MPI_INC} $(CUDA_INCL) ${ADIOS_INC} -I${FFTW_INC}

#[ Define library paths.
LFLAGS = -L../lib -L${MPI_LIB} $(CUDA_LIBS) ${ADIOS_LIB} -L${FFTW_LIB} $(FFTW_LIB_FLAGS)

#[ Define any libraries to link into executable:
#[   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#[   option, something like (this will link in libmylib.so and libm.so:
LIBS = 

#[ Define the C source files.
SRCS = alloc.c ffts.c initialization.c io_adios.c io_tools.c comms.c mugy.c utilities.c data.c finalize.c
ifeq ($(USE_GPU),true)
    SRCS_DEV = initialization_dev.cu alloc_dev.cu data_dev.cu
else
    SRCS_DEV =
endif

#[ Define the C object files.
#[ Uses Suffix Replacement within a macro:
#[   $(name:string1=string2)
#[         For each word in 'name' replace 'string1' with 'string2'
#[ Below we are replacing the suffix .c of all words in the macro SRCS
#[ with the .o suffix
OBJS = $(SRCS:.c=.o)
ifeq ($(USE_GPU),true)
    OBJS_DEV = $(SRCS_DEV:.cu=.o)
else
    OBJS_DEV =
endif

#[ Define the executable file.
MAIN = mugy.x
ifeq ($(USE_GPU),true)
    MAIN_DEV = mugy_dev.o
else
    MAIN_DEV =
endif

#[ The following part of the makefile is generic; it can be used to
#[ build any executable just by changing the definitions above and by
#[ deleting dependencies appended to the file from 'make depend'

.PHONY: depend clean

all:    $(MAIN)
	@echo  Source mugy.c has been compiled into executable: mugy.x

$(MAIN): $(MAIN_DEV) $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(MAIN_DEV) $(OBJS) $(OBJS_DEV) $(CUDART) $(LFLAGS) $(LIBS) -o $(MAIN)

#[ This is a suffix replacement rule for building .o's from .c's
#[ it uses automatic variables $<: the name of the prerequisite of
#[ the rule(a .c file) and $@: the name of the target of the rule (a .o file)
#[ (see the gnu make manual section about automatic variables).
$(MAIN_DEV): $(OBJS_DEV)
	$(NVCC) $(CUDAFLAGS) -dlink $(INCLUDES) $^ -o $(MAIN_DEV)

$(OBJS_DEV): %.o: %.cu
	$(NVCC) $(CUDAFLAGS) $(INCLUDES) -dc $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

#[ DO NOT DELETE THIS LINE -- make depend needs it.
