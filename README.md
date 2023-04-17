# mugy
Multiscale gyrofluid code

Written in C and CUDA, with postprocessing scripts in Python.

Installation
------------

The code has the following dependencies:
- GSL (GNU Scientific Library)
- CUDA
- MPI
- FFTW
- ADIOS2

Follow these steps to build mugy:
1. Specify paths to include and library directories, as well as some flags, for each of these dependencies in the ```setup.sh``` script.
2. Run the setup script: ```. ./setup.sh```.
3. Run the makefile: ```make```.

This produces the executable ```mugy.x```.


Runing mugy
-----------

mugy must be run with MPI and requires at least 2 command line inputs, the full path to the input file and the full path to the output directory. For example, to use the input file ```mugy_test.in``` and output to ```/scratch/jdoe/test/``` use:

```mpirun -np 1 ./mugy.x ./mugy_test.in /scratch/jdoe/test/```

An optional third command line argument may be specified: the restart directory (the directory of a previous simulation that we wish to continue). For example, if the restart directory is ```/scratch/jdoe/start/``` then this restart simulation may be launched with

```mpirun -np 1 ./mugy.x ./mugy_test.in /scratch/jdoe/test/ /scratch/jdoe/start/```


Developer notes
---------------

Please try to use following naming conventions whenever possible: 
- Datatype names must end with ```_t```.
- ```enum``` options must be capitalized.
- append ```mugy_``` as a prefix to structs, datatypes, and functions that may be publicly visible.
