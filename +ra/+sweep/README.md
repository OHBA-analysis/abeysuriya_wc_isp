### Using the sweep system

Sweeps are implemented using the `Pepper` repository. To run them locally, do the following

First, make the simulation output file

	p = ra.sweep.simulate(ra.network.line(2),ra.unit.deco2,linspace(1,80,3),linspace(0,0.5,3));
	o = p.run();
	save('temp','-struct','o')

Note that `temp` has the same structure as the results file that would be generated if the simulation had been performed on the cluster. Next, process the file using `analyse`

	p = ra.sweep.analyse('temp')
	o = p.run();
	save('temp_results','-struct','o')

Lastly, display the results file

	ra.sweep.display('temp_results')

### General overview

The sweep is implemented using the `ra.sweep.simulate` which inherits from `Pepper`. Analysis is delegated to a secondary run using `ra.sweep.analyse` which enables the analysis to be re-run in parallel without having to run the simulation again. The implementation is as follows

- `simulate` stores member variables for the unit, network, and coupling/velocity values to test
- `simulate.get_inputs` returns matrices of coupling/velocity
- `simulate.do_work` calls `ra.wc_metastable.template` to run the simulation and return a `WilsonCowan` object
- `simulate.postprocess` reshapes the `WilsonCowan` array into a 2D grid, the same as the velocity and coupling matrices 
- `simulate.assemble` first saves the output to `results.mat` by calling `Pepper.assemble()`, but then performs a system call to the `submit` script to immediately issue an analysis job for the same run
- `analyse` stores the file name of the input result file to be read in. The working directory is still specified by the `submit` script
- `analyse.get_inputs` loads the original result file
- `analyse.do_work` returns a struct with the computed measures for a single run
- `analyse.postprocess` reshapes the output arrays into matrices
- `analyse.assemble` first saves the output to `results.mat` by calling `Pepper.assemble()`, and then moves the file to the directory for the original job 