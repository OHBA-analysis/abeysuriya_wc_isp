## README

Code to perform simulations in 

**Abeysuriya et. al. "A biophysical model of dynamic balancing of excitation and inhibition in fast oscillatory large-scale networks"**

### Project structure

The project is built on two levels

- `integrate_wc.m` is an integrator that takes in all the parameters of the model/system and integration settings, and outputs the solution of the system. To perform basic simulations or to build on the work, we would suggest starting here.
- The code in the `+ra` folder which facilitates constructing the inputs for `integrate_wc`, scanning over parameter values, and performing the analysis. Note that this code was built to integrate into a larger simulation engine, and it is thus somewhat more complicated. We envisage the primary use of this code as being to perform the simulations shown in the paper, or to make minor changes (e.g. different parameter values). 

### Setup

#### Prerequisites

Basic simulation 
- Matlab (tested on R2016a and R2017a). 

Whole-brain simulations and analysis
- Matlab Signal Processing Toolbox and Stats toolbox
- Fieldtrip - To set up Fieldtrip, decompress `data_files/fieldtrip.tar.bz2`, add that folder to the Matlab path, and run `ft_defaults` 
- MEG-ROI-nets - To set up ROInets, decompress `data_files/roinets.tar.bz2` and add the `MEG-ROI-nets` folder to the path

You could alternatively satisfy the Fieldtrip and MEG-ROI-nets dependencies by installing OSL (see [https://ohba-analysis.github.io/osl-docs/](https://github.com/OHBA-analysis/) or [https://github.com/OHBA-analysis/](https://github.com/OHBA-analysis/)), which also provides a range of functionality to facilitate plotting and analysis of experimental data. 

#### Startup

Add the main folder (the one containing `+ra`) to the path, or just work in that directory.

### Usage

See `examples.m`. Note that the ordering of the ROIs in the original code and in the paper are different. Image plots
need to use `ra.analysis.reorder_matrix` to change the ordering. The supplementary material online has had
this reordering already applied. The files provided with this code (in the `data_files` folder) use the original ordering. 

The original simulations were performed in Matlab with a MEX-compiled C++ simulation engine. The implementation here is a more recent pure Matlab implementation that we have verified produces identical results to within numerical precision. This engine scales differently to the C++ implementation, and is significantly faster than the original code for whole-brain networks. Because the most interesting parameters result in chaotic dynamics, these numerical differences cause the timecourses to come out differently after a moderate period of time (on the order of a few seconds) because transitions between oscillatory states occur at slightly different times, which results in the states diverging rapidly after the first time this occurs. The same effect occurs if the original simulation code is compiled and run on a different platform, although it takes a bit longer (tens of seconds) for the difference to become noticable. While this does affect the ability to produce the exact output as used in the paper, the situation is the same as if a different random seed were used for the noise. The saved simulations for the results in the paper are available on request (they are ~200GB in total and have thus not been included here).

Published simulations sweeping over delay and coupling values were performed using the code in `+ra/+sweep` as called by the `sweeps` folder. These map as follows

- `deco_rk4.m` - Simulations without ISP e.g. Fig 4a
- `deco_rk4_adjustmean` - Simulations with coarse balancing e.g. Fig 9a
- `deco_rk4_isp` - Simulations with ISP e.g. Fig 4b
- `deco_rk4_isp_noise_test` - Simulations with different noise levels e.g. Fig. S2 in Supplementary Material
- `deco_rk4_isp_target_10_hd` - Simulations with low ISP target e.g. Fig. 8a
- `deco_rk4_isp_target_30_hd` - Simulations with high ISP target e.g. Fig. 8b

These files are included for reference, but they were called by additional infrastructure specific to the cluster that was used, so won't run directly. See `examples.m` for how to perform the type of simulations used in the paper using the code provided here. 

### Included files

- `cortical_network.mat` stores the network connectivity. This is the same as the separate files in the supplementary online material.
Note that the ordering is different to in the paper. The ROI names are stored internally, e.g.

    net = ra.network.import('cortical_network.mat')
    net.nodes

Note also that the `Coordinates` correspond to 8mm MNI coordinates and match those in the parcellation file 
included in the supplementary material, notwithstanding the reordering.
- `dk_cortical.nii.gz` - NIFTI file of the parcellation, with ordering that matches `cortical_network.mat`
- `dk_cortical.txt` - ROI labels, with ordering that matches `cortical_network.mat`
- `data_fc.mat` - The subject-specific alpha band connectivity matrices from data
- `fig_5_simulation` - Precomputed results for the simulation shown in Fig. 5
- `isp_means.mat` stores the mean value of `cie` used for the coarsely balanced simulations in Fig. 9a.

### Key objects

There are three key classes that are used by this code, detailed below:

#### TimeSeries

The `TimeSeries` object provides a range of useful methods for working with time series. See the code for full details, but key functionality includes

- `resample` to change the sampling rate
- `trim` to remove time from the start or end
- `filter` to apply a filter
- `orthogonalize` to run the signal orthogonalization routine
- `envelope` to compute amplitude envelope and phase timecourses. Can also accept filter and orthogonalization arguments 

See also `examples.m` for how to use this functionality as part of the analysis. 

#### Model

The `Model` object provides a standard interface to run simulations of networked oscillators, although only the Wilson-Cowan model is provided here. The key usage is

- `wc = ra.model.WilsonCowan(u,net)` creates a model object, where `u` is a unit struct (e.g. returned by `ra.unit.deco_10`) and `net` is a network object
- `wc.units` stores the unit parameters
- `wc.network` stores the network object
- `wc.params` stores network level parameters
- `wc.options` sets some options for the integration, including integration time (`wc.options.tspan`) and time step (`wc.options.tstep`) 
- `wc.run()` runs the simulation. `wc.run('burn',[t1 t2],'downsample',fs)` will downsample the output to `fs` (including an antialiasing filter) and then trim `t1` seconds from the start, `t2` seconds from the end. Best used together, to remove filter edge artefacts.
- `wc.result` contains the output for all nodes
- `wc.excitatory` and `wc.inhibitory` contain the outputs for just the excitatory and inhibitory nodes, respectively
- `wc.extend()` initializes the simulation with the current output, thus extending the simulation run
- `wc.plot()` will plot the excitatory timeseries, for convenience

#### Network

The `Network` object creates different types of networks - graphs where nodes correspond to brain regions. Each edge has a weight and a distance. Each node has a name and associated coordinates (which are of course only meaningful in the case of real brain networks). 

- `net = ra.network(cmat,dmat)` creates a network given a square connection matrix and distance matrix
- `net = ra.network.null(n)` creates a network with `n` nodes and no edges
- `net = ra.network.uniform(n)` creates a fully connected network with `n` nodes and all possible connections between them (no self connections)
- `net = ra.network.import(n)` loads a network from a file (specifically, we provide `data_files/cortical_network.mat` which is the network used in the paper) 
- `[cmat,dmat] = net.netmats()` returns connection and distance matrices. Note that the network object has fields `net.edges` which is an adjacency-list representation of the graph, and the `netmats()` function converts it to a matrix representation. 


### Folder structure

- `+ra` 
	- `+analysis` Older analysis scripts
	- `+data` - Functions for processing experimental data
	- `+model` - Objects for constructing and running simulations
	- `+paper` - Functions used to generate figures in the paper
	- `+plot` - Various plotting routines
	- `+sweep` - Functions used to perform long runs and analysis for the paper
	- `+unit` - Operating point parameters for Wilson-Cowan
	- `+utils` - Helper scripts
	- `+wc` - Some helper functions for working with Wilson-Cowan model objects
	- `@network` - Implementation of weighted network container class
- `sweeps` - Contains driver scripts used to perform runs for the paper
- `data_files` - Additional content, described above in 'Included files'

