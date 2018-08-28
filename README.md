### N-Body evolution via 4th-order Hermite Integration ###
----------------------------------------------------

This FORTRAN 90 repository sets up N gravitating bodies and solves their equation
of motion numerically via the 4th-order Hermite algorithm, with a global (adaptive) timestep.  

The code also calculates orbital elements (by default around the centre of mass).  The principal force is
Newtonian Gravitation - no other forces are implemented.

### Outputs ###
---------------

The code can produce two types of output:

1) One file per particle, to easily plot particle properties vs time

2) Snapshot files, where all particle data at a single timestep is recorded


### Compilation and Execution ###
----------------------------------

This code was developed and tested on gfortran 8.2.0, and compiled with standard Makefile software.

The code is compiled using `src/Makefile` contained within the repo, i.e. to compile
`> cd src/`

`> make`

Which produces the `nbody_hermite` executable.  The code is then run with the command

`> ./nbody_hermite <paramfile> `

Where <paramfile> is the parameters file needed to setup the N Body system - for an example see `src/nbody_hermite.params`
  
### Plotting ###
----------------

The `plot/` subdirectory contains a set of Python scripts (mostly Python 3.0) to plot both single files, snapshots and the log files that are outputted by the code.  Most of the heavy lifting (data reads, setting up figures) is done in the `io_nbody_hermite` module, which gives useful functions for writing more sophisticated plotting scripts.
