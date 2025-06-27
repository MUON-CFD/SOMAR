Welcome to the SOMAR repository!
=====

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14609305.svg)](https://doi.org/10.5281/zenodo.14609305)

SOMAR stands for The Stratified Ocean Model with Adaptive Refinement. It is free software using the [LGPL license](https://www.gnu.org/licenses/lgpl-2.1.html 'The GNU Lesser General Public License, version 2.1 applies.') and provided jointly by [Thomas Jefferson University's College of Humanities and Sciences](https://www.jefferson.edu/academics/colleges-schools-institutes/humanities-sciences.html 'TJU\'s CHS website') and [Arizona State University's School of Engineering of Matter, Transport, and Energy](https://semte.engineering.asu.edu/ 'ASU Engineering website').


Main Features
-------------
- **Nonhydrostatic made fast -** Blending [the leptic method](https://doi.org/10.1016/j.jcp.2011.06.022) with semicoarsened multigrid, SOMAR solves the Boussinesq Navier-Stokes equations *without* the hydrostatic approximation in order to properly model the internal waves and tides that are ubiquitous in the ocean.


- **Anisotropic adaptive mesh refinement (AMR) -** A coarse underlying grid along with on-the-fly local refinement of transient features eliminates unnecessary computation in most of the domain. Refinement can occur in some or all directions by varying amounts. Also, refinement occurs in both space *and* time to provide a drastic speedup of computation and reduction of memory usage.


- **Separation of background density and its deviation -** By splitting the density field into a vertical background stratification and a deviation, we relieve the Poisson solver of computing the associated hydrostatic component of the pressure. This treatment, already implemented in some regional models including [MITgcm](http://mitgcm.org/ 'The MITgcm website'), also prevents diffusion of oceanic features that are maintained by unmodeled phenomena.


- **Large Eddy Simulation (LES) -** At local turbulent hotspots, AMR helps resolve the bulk of the energy cascade and LES parameterizes the rest. This two-way feedback between models uses the subgrid closure scheme of [Ducrose, et. al.](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/largeeddy-simulation-of-transition-to-turbulence-in-a-boundary-layer-developing-spatially-over-a-flat-plate/C277DE968A1FD929D3CB05FDBC434AAD)


- **Familiar methods -** Unlike SOMAR v1.0, this newer version uses methods familiar to physical oceanographers such as Arakawa-C grids and Runge-Kutta time integration. We also avoid upwinding by using simple, centered finite differences for the advection terms.


- **Error control -** Embedded RK schemes are supported. Users can limit the velocity error using a variety of methods (local, PI, or PID controllers). From the user's perspective, all that is needed is an error tolerance specified in a text-based input file and SOMAR will take care of the rest.


- **Python post-processing -** Analyzing SOMAR's output can be accomplished in [VisIt](https://visit-dav.github.io/visit-website/index.html 'The VisIt website'), [ParaView](https://www.paraview.org/ 'The ParaView website'), or with simple Python scripts. Python can also be used to create initial conditions, create custom force functions, and process data on-the-fly.


- **Highly parallelizable -** SOMAR uses MPI and is build on the [Chombo framework.](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)


Work is currently underway to incorporate the immersed boundary method into SOMAR. This will allow us to incorporate complex boundaries provided by digital elevation models. This effort is largely complete and will become available once its associated manuscript is published.


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
Up-to-date documentation is available at [https://somar.readthedocs.io](https://somar.readthedocs.io).
