Welcome to the SOMAR repository!
=====

SOMAR stands for The Stratified Ocean Model with Adaptive Refinement. It is free software using the [LGPL license](https://www.gnu.org/licenses/lgpl-2.1.html 'The GNU Lesser General Public License, version 2.1 applies.') and provided jointly by [Thomas Jefferson University's College of Humanities and Sciences](https://www.jefferson.edu/academics/colleges-schools-institutes/humanities-sciences.html 'TJU\'s CHS website') and [Arizona State University's School of Engineering of Matter, Transport, and Energy](https://semte.engineering.asu.edu/ 'ASU Engineering website').

Release info
-----
The latest stable versions have been tagged as:
- **Version 1.2-alpha1**: Introduces static mesh stretching and works with any number of levels.
- **Version 1.1-alpha1**: Tested on simple Cartesian grids and 2 levels.

Work is currently underway to incorporate the immersed boundary method of Charles Peskin into SOMAR to incorporate complex boundaries.


Features
-----
- **Nonhydrostatic -** Our model solves the Boussinesq Navier-Stokes equations *without* the hydrostatic approximation in order to properly model the internal waves and tides that are ubiquitous in the ocean.


- **Separation of background density and its deviation -** By splitting the density field into a vertical background stratification and a deviation, we relieve the Poisson solver of computing the associated hydrostatic component of the pressure. This treatment, already implemented in some regional models including [MITgcm](http://mitgcm.org/ 'The MITgcm website'), also prevents diffusion of oceanic features that are maintained by unmodeled phenomena.


- **Anisotropic grid refinement -** A coarse underlying grid along with dynamic local refinement over transient features eliminates unnecessary computation in large portions of the domain. Since the background stratification requires some level of vertical resolution, it is often the case that we only need further resolution in the horizontal. Our anisotropic refinement methods are capable of providing additional cells only in those directions that are under-resolved. Furthermore, our coarse grids operate on larger timesteps than the finer grids. This refinement in both time and space provides a drastic speedup of computation, minimizes the number of required Poisson solves, and ensures all levels are evolving at a Courant number close to one, reducing numerical dissipation.


- **Fast multigrid Poisson solver -** In the absence of the hydrostatic approximation, we are faced with an exceedingly expensive Poisson problem for the pressure. To relieve the computational burden, SOMAR is equipped with a multigrid solver that can reduce the condition number of anisotropic problems using semi-coarsening. The most effective coarsening strategy is created automatically, without asking the user to tune parameters.


- **Large Eddy Simulation -** AMR helps refine local hotspots, but we often do not wish to resolve the flow down to the dissipative scales. For this, we turn to the closure scheme described by [Ducrose, et. al.](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/largeeddy-simulation-of-transition-to-turbulence-in-a-boundary-layer-developing-spatially-over-a-flat-plate/C277DE968A1FD929D3CB05FDBC434AAD)

- **Familiar methods -** Unlike SOMAR v1.0, the newer versions use methods familiar to physical oceanographers such as Arakawa-C grids and Runge-Kutta time integration. We also use simple, centered finite differences for the advection terms rather than upwinding.

- **Highly parallelizable -** SOMAR uses MPI and is build on the [Chombo framework.](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)


<!-- Software prerequisites
-----
TODO


Compilation
-----
TODO


Taking a test drive
-----
Before creating your own simulation, you should successfully test one of the packaged simulations that come with the solver. In the `exec` folder, you will see several input files. Each of these files are targeted to run a specific problem on a specific machine and are identified accordingly (`inputs.problem.machine`). As a rule, an input file for machine A should never be altered to work on machine B. So if we want to run, say, the lock exchange demo problem on a machine called HAL, we should first copy one of the lock exchange input files to `inputs.LockExchange.HAL`, and then edit it as we please.

Open this new file and look for the following lines (it's okay if some of them are missing).

```
# amr.restart_file = chkpt_000010.2d.hdf5
# plot.plot_prefix = plot_
# plot.checkpoint_prefix = chkpt_
# plot.plot_period = 0.1
plot.plot_interval = 1
plot.checkpoint_interval = 100
```
In this listing, several of the lines are commented out with a `#`. This usually means we either do not want a particular feature or we are happy with the default values. In this case, we do not want to restart a simulation from a saved state and we are happy with the default file name prefixes. As you can see, all plot files will be prefixed with `plot_` and all checkpoint files (files used to restart a simulation) will be prefixed with `chkpt_`. If we want the plot files to go into a different directory, we can uncomment the plot\_prefix line to read `plot.plot_prefix = /home/user/myPlotFolder/plot_`. Other input parameters can similarly be altered as you see fit.

With your code compiled and input file prepared, you are now ready to run the lock exchange demo. To run the demo in serial, use `./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. To run this simulation in parallel over 8 processors, use `mpirun -np 8 ./somar2d.[config].OPTHIGH.MPI.ex inputs.LockExchange.HAL`. The result should be two sets of [HDF5](http://www.hdfgroup.org/HDF5/ 'The HDF group website') files -- one set of checkpoint files that are used to restart a simulation, and one set of plot files that can be viewed in [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit 'The VisIt webpage'). -->


Documentation
-----
Documentation will eventually be hosted on [ReadTheDocs.org](https://readthedocs.org/), but for now you can compile the docs locally. First, make sure [sphinx](https://www.sphinx-doc.org/en/master/) is installed on your machine. You can typically do this via `pip3 install sphinx` or `apt-get install python3-sphinx`. Once installed, switch to SOMAR's `docs/sphinx` folder and run `make html`. View the docs by pointing your web browser to `docs/sphinx/build/html/index.html`.
