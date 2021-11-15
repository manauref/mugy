#[
#[ Makefile for mugy
#[
#[
#[ Users need to set the following environment variables:
#[
#[      CCOMPILER, ADIOS_DIR
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



#[ Define the C compiler to use.
CC = /Users/manaure/Documents/multiscale/code/openmpi4.1.1-intel19.1/bin/mpicc

#[ ADIOS flags:
ADIOS_FLAGS = -DADIOS2_USE_MPI
ADIOS_LIB = -Wl,-rpath,/Users/manaure/Documents/multiscale/code/adios2-openmpi4.1.1-intel19.1/lib /Users/manaure/Documents/multiscale/code/adios2-openmpi4.1.1-intel19.1/lib/libadios2_c_mpi.2.7.1.dylib /Users/manaure/Documents/multiscale/code/adios2-openmpi4.1.1-intel19.1/lib/libadios2_c.2.7.1.dylib
ADIOS_INCL = -I/Users/manaure/Documents/multiscale/code/adios2-openmpi4.1.1-intel19.1/include

#[ Define any compile-time flags.
CFLAGS = -Wall -std=c99 $(ADIOS_FLAGS)

#[ Define any directories containing header files.
INCLUDES = -I../include $(ADIOS_INCL)

#[ Define library paths.
LFLAGS = -L../lib $(ADIOS_LIB)

#[ Define any libraries to link into executable:
#[   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#[   option, something like (this will link in libmylib.so and libm.so:
LIBS = 

#[ Define the C source files.
SRCS = alloc_mugy.c initialization.c io_tools.c mpi_tools.c mugy.c utilities.c data_mugy.c 

#[ Define the C object files.
#[ Uses Suffix Replacement within a macro:
#[   $(name:string1=string2)
#[         For each word in 'name' replace 'string1' with 'string2'
#[ Below we are replacing the suffix .c of all words in the macro SRCS
#[ with the .o suffix
OBJS = $(SRCS:.c=.o)

#[ Define the executable file.
MAIN = mugy.x

#[ The following part of the makefile is generic; it can be used to
#[ build any executable just by changing the definitions above and by
#[ deleting dependencies appended to the file from 'make depend'

.PHONY: depend clean

all:    $(MAIN)
	@echo  Source mugy.c has been compiled into executable: mugy.x

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

#[ This is a suffix replacement rule for building .o's from .c's
#[ it uses automatic variables $<: the name of the prerequisite of
#[ the rule(a .c file) and $@: the name of the target of the rule (a .o file)
#[ (see the gnu make manual section about automatic variables).
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

#[ DO NOT DELETE THIS LINE -- make depend needs it.
