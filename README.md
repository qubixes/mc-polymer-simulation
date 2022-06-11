# Monte Carlo Polymer simulation


## General note:

This repository contains code used in several [papers](https://www.sciencedirect.com/science/article/pii/S0021999118301098) for simulating polymers using
Monte Carlo simulation. The repository contains a few different implementations, that are interoperable through the polymer storage files.

As a general note: I am not working in polymer physics anymore, so the code hasn't been maintained in the last three years (as of 2022).

It consists mainly of C/cuda/OpenCL executables and Bash scripts that tie it all together. Generally, it is not necessary to directly run the executables (some of which have 10+ arguments in a specific order). Instead, run the scripts. 

To run the package there are the following prerequisits (some are optional):

- *nix (Linux and OS X should work fine)
- C/cuda/OpenCL compiler
- gnuplot [visualization, aqua terminal is usually used (change to wxt for Linux)]
- povray [3D visualization]
- Fortran 90 compiler [secondary structure analysis]
- Bash

Unfortunately, there are no real automake/install scripts, so you might need to add some directories to Makefiles if your compiler cannot automatically find the necessary libraries.



## Directories

### analysis

This directory contains programs to analyze data.

### cpupol

This directory contains ELP (Elastic Lattice Polymer model) programs. They are not particularly efficient, but relatively easy to program.

### denspol

The DLP (Dense elastic Lattice Polymer model) programs are stored there. Another CPU implementation, but in this case it contains a novel way to create very high density melts, effectively up to 5 times higher than the reference
implementation. This gives very big benefits in simulating melts in the reptation regime, where it even outperforms the GPU implementations
(1 CPU core vs 1 GPU). Unfortunately, the paper has ended as a (mostly complete draft) and was never published.

### gpupol

It contains the old GPU program for both linear and ring polymers. For ring polymers this is pretty much legacy code, however it contains the only GPU program for linear polymers.

### gpupol2

The new GPU programs are stored here (versions 2 and 3). In general, the last version should be the fastest of the three.

### scripts

This directory contains all sorts of scripts, mostly visualization ones. 

### secstruct

Secondary structure analysis with Daniel Jost's program. It contains a wrapper to apply it to the TUV data format. 

### util

Some C files that are used with different programs such as RNG/IO.



### Scripts

### do_run.sh

This is the main script to run simulations. While you can run it without any arguments for testing, there are many arguments to define the simulation parameters. They are not all guarenteed to work in each and every combination. Each program might have different features, which are simple not implemented on others (or near impossible to implement). For the current default parameters, that can be checked directly in the file. The important arguments are as follows (compatibility is signified with parentheses):

`-x|--exec [base executable name]`

	Examples are "efpol", "gpupol", "denspol". Different variants such as OpenCL/Cuda or linear/ring polymers have the same base executable name.
	
`-r|--basedir [directory]`

	Base directory where the simulations are found.
	
`--outdir [directory]`

	Directory where the current simulation is to be stored. Overrides -r argument.

`-b|--bend [bending energy] (denspol)`

	Set the bending energy between two consecutive bonds. 
	
`-u|--double [number of generations]`

	For hierarchical building of systems (with relaxation). If number of generations equals zero, no hierarchical building is done.

`-g|--grid [lattice size]`

	Length of the edge of the simulation box in number of lattice sites. Total number of lattice sites is this cubed. Not every size is permisible for the GPU algorithms.

`-s|--seed [seed]`

	Seed of the simulation. Integer value <2^32.

`-t|--time [time]`

	Number of time steps of the simulation.

`-n|--nmono [nmono]`

	Number of monomers of each polymer in the simulation. Generally, the number of polymers is adjusted to density and lattice size.

`-i|--interval [interval]`

	The number of time steps between writing the results to a file. Generally, it is expected that [time] is a multiple of [interval]. Most of the time writing between 300 and 3000 files results in a time series that is reasonable to analyze.

`--ring`

	Simulate ring polymers.

`--linear`

	Simulate linear polymers.

`-d|--density [density]`
	Density of the system in number of monomers per lattice site. Normal ELP simulations use about 1.2, while DLP is normally around 7.2.

`-f|--fast-eq (-)`
	Try to equilibrate faster by redistributing stored length. Never worked right/obsolete. 

`-m|--multilength [file] (denspol)`

	To enable polymers with different lengths supply a file with different lengths and topology. If the harmonic potential is enabled (see further on), then the file should also contain the contacts/pairs of monomers that attract each other. This is mainly used for reconstruction.

`--short (all?)`

	A way to produce equilibrated simulations with shorter time intervals. Needs the long equilibrated simulations first as a starting configuration (detection is automatic).

`--topo|--chain-crossing (denspol)`

	Allow chain crossing in the simulation. Statics should remain the same (dynamics obviously not!)

`--lattice-sphere (denspol)`

	Instead of periodic boundary conditions with no obstacles, start with a confined sphere. Actual radius is a function of the lattice size (not  exactly at the boundaries, check what suits you).

`--hp|--hp-strength (denspol)`

	Strength of the harmonic potential between monomers. Is used for reconstruction.


### check_similarity.sh

This file compares a reconstruction to the original. Arguments (all mandatory):

- `original directory`

	This can be a DLP/ELP directory, or one with a configuration supplied by Pascal. In the ELP format, it will compare the last file (which is the one supposedly reconstructed). 

- `reconstruction directory`

	A DLP/ELP directory that contains the reconstruction.

- `boundary conditions`

	To enforce no periodic boundary conditions, supply "static" (without the quotes).

It outputs directly to stdout the error as a function of time. In case of DLP data, the format is:

`t err_int_forw err_ext_forw err_int_back err_ext_back err_contact`

"int" stands for internal distances within polymers, while "ext" is for distances between different polymers. The contact error is just the distance between two monomers that are in contact. Forward or backward errors should not be much of a difference. The difference between the two is which of the original or reconstructed configurations is interpolated.


### create_recon_batch.sh, psmn_reconstr.sh

These scripts together create batches for reconstruction purposes, with different parameter: number of contacts, chain-crossing (or not), potential strength. Each of the parameter sets is run with 8 different seeds. 

- `[queue name]`
	It takes the name of the queue as its argument.



### desec.sh, exec.sh, psmn_denspol_ring.sh

These basically contain examples of how to use the scripts.



### reconstruct.sh

This script does the reconstruction of systems, either pascal data or DLP data. In principle it can be done solely with the do_run.sh script, but it does a lot of bookkeeping that makes life a lot easier. There are less options to choose from than in do_run.sh (some are just not necessary). There is one new option:

`--nsamples [number of contacts]`

	This allows for setting the number of contacts that are taken account for the reconstruction. If it is less than the actual number of contacts, random discarding is used.

The output is in a new directory `[old_dir]/../reconstruction/[directory depending on parameters]_b[generation]`.


### update.sh

This script runs the analyses. It only runs the ones that do not have their results. If you want to force a rerun of the analysis of a certain directory, either "touch simulation_settings.txt" or delete the result file that you want to update (it is faster to analyze individual results). 


### Makefile

This is mostly the standard stuff:

- all

	Compile all executables.

- install 

	Compile all executables and install them in the ./bin directory.

- clean

	Remove all executables/object files/some temporary files

- tar

	First clean everything, then make tar.gz file and put it in ../tar. It ignores any ./data directory, so it does not get included in the tar.gz file.

The "RELEASE" macro is used throughout to keep track of which version is used for simulation. This helps in case of a bug to see if it is necessary to rerun.


- git

There is a git tree embedded in the tar files. 
