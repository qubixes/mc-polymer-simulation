# Description of scripts

This file describes and documents the scripts in this directory. The vast majority of these scripts have to do with plotting or otherwise visualizing the results (run update.sh in case of missing files). 

## util.sh

This file is included in many of the other files, and is basically used as a library. The most important functions are:

- `get_attr [attribute] [directory]`

	The attribute here is one in the file "simulation_settings.txt", which is always supposed to be there.

- `get_title [dir] [format string]`

	Replaces "%n" with the lenght of the polymers in a string.

- `needs_update [file_orig] [file_dest]`

	If `[file_dest]` does not exist or is older than `[file_orig]`, it returns 0, otherwise 1.

- `cross_sect_files [fileA] [fileB]`

	Reads two files and outputs one file with the first column in common.

- `get_last_tfile [dir]`

	Finds the name of the last tuv file in a directory (directory is not included in the name).

- `get_last_savfile [dir]`

	Same but for the lattice file for the denspol algorithm.

- `get_last_t [dir]`

	Get the highest timestep value of file in the directory.

- `get_dirs [args]`
	This function has similar arguments as the do_run.sh script (but lacking many options). The arguments are a filter for the directories returned. E.g. with --linear, only simulations with linear are returned. For some of the arguments it can be submitted multiple times, resulting in both being accepted. Giving -n 100 -n 200, should return all directories with either 100 or 200 monomers.
	
	`-n|--nmono [nmono] (multiple)`
	
		Number of monomers.
	
	`-d|--density [density] (multiple)`
	
		Density of the system.
	
	`-b|--double [steps] (multiple)`
	
		Number of hierarchical building steps.
	
	`-l|--linear`
	
		Only linear polymers.
	
	`-r|--ring
	
		Only ring polymers.
	
	`-e|--equilibrated`
	
		Only equilibrated systems.
	
	`-x|--exec [algorithm] (multiple)`
	
		The algorithm used for the system, e.g. `gpupol[,2,3]/denspol/efpol`.
	
	The algorithm sorts the directories with increasing length and outputs them to standard output.


## copy_dbl.sh, double_exec.sh fix_attr.sh get_dead_dirs.sh

These are just small scripts to do some specific running copying of things. Can be deleted for other projects. Not recommended to run without knowing what it does.


## get_dirs.sh

A version of get_dirs in util.sh. Can be useful in C programs. Otherwise, just use the one in util.sh.


## list_runs.sh

This script lists the directories 

## Makefile

For removing old temporary files, and also it contains a shortcut to create .eps files from .tex gnuplot output files.

## plot_cmsdif.sh

Plot the center of mass displacement as a function of time. Takes [get_dirs] arguments (see util.sh). 

## plot_corcms.sh (legacy?)

Plot the correlation of cms velocity.

## plot_diff.sh

Plot the diffusion coefficient as a function of z=N/Ne. [get_dirs] type arguments. Second plot is one that gives a measure of equilibration time needed for the systems. This measure is based on diffusing the center of mass over its own radius of gyration (times a prefactor). You have to manually adjust the type of polymer to get the right control lines in the plot. 

## plot_dist_rec.sh

Plot a detailed comparison between Pascal's data and its reconstruction. The first plot shows the average and standard deviation of the distance in Pascal's data as a function of the distance in the reconstruction. The second plot is the same, but with normalized distances. The third plot shows the monomer density profile as a function of the distance from the center of mass. Readjust the normalization for maximum visibility.

## plot_double.sh

This is a program to plot the radius of gyration of polymer size after relaxing for $T time. The value for $T is set directly in the program. It takes [get_dirs] arguments. 

## plot_efficiency.sh

This script takes no arguments, but reads from the file "ne_list.dat" instead. The first plot shows the efficiency in terms of how many Monte Carlo steps are needed to reach $\tau_e$. The second plot is corrected for how many Monte Carlo Steps can be done per second with each algorithm. On the x-axis is the density of the systems. It is just there so that they are easier to distinguish. There are differences in models with different bending energies/stored length energy, which is why the curves are not smooth.

## plot_eq.sh

This plots the equilibration of systems using the radius of gyration as a function of time. The second plot shows the radius of gyration as a function of time. The second plot shows the convergence to the average as a function of the average center of mass displacement. 

## plot_genom.sh

It plots the overlap parameter as a function of N, the genomic distance as a function of N, and an estimate of Ne (times a constant). The estimate is made by finding the crossing point of the |i-j|, |i-j|^{2/3} lines. Input as usual in [get_dirs] form. 

## plot_magdip.sh

Plotting the magnetic dipole moment/magnetic radius in varius ways. Input as in [get_dirs]. 
