Welcome to the SOMAR repository!
=====

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14609305.svg)](https://doi.org/10.5281/zenodo.14609305)

SOMAR stands for The Stratified Ocean Model with Adaptive Refinement. It is free software using the [LGPL license](https://www.gnu.org/licenses/lgpl-2.1.html 'The GNU Lesser General Public License, version 2.1 applies.') and provided jointly by [Thomas Jefferson University's College of Humanities and Sciences](https://www.jefferson.edu/academics/colleges-schools-institutes/humanities-sciences.html 'TJU\'s CHS website') and [Arizona State University's School of Engineering of Matter, Transport, and Energy](https://semte.engineering.asu.edu/ 'ASU Engineering website').


Main Features
-------------
- **Nonhydrostatic made fast -** Blending [the leptic method](https://doi.org/10.1016/j.jcp.2011.06.022) with semicoarsened multigrid, SOMAR solves the Boussinesq Navier-Stokes equations *without* the hydrostatic approximation in order to properly model the internal waves and tides that are ubiquitous in the ocean.

- **Anisotropic adaptive mesh refinement (AMR) -** A coarse underlying grid along with on-the-fly local refinement of transient features eliminates unnecessary computation in most of the domain. Refinement can occur in some or all directions by varying amounts. Also, refinement occurs in both space *and* time to provide a drastic speedup of computation and reduction of memory usage.

- **Separation of background density and its deviation -** By splitting the density field into a vertical background stratification and a deviation, we relieve the Poisson solver of computing the associated hydrostatic component of the pressure. This treatment, already implemented in some regional models including [MITgcm](http://mitgcm.org/ 'The MITgcm website'), also prevents diffusion of oceanic features that are maintained by unmodeled phenomena.

- **Large Eddy Simulation (LES) -** At local turbulent hotspots, AMR helps resolve the bulk of the energy cascade and LES parameterizes the rest. This two-way feedback between models uses the subgrid closure scheme of [Ducros, et. al.](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/largeeddy-simulation-of-transition-to-turbulence-in-a-boundary-layer-developing-spatially-over-a-flat-plate/C277DE968A1FD929D3CB05FDBC434AAD)

- **Familiar methods -** Unlike SOMAR v1.0, this newer version uses methods familiar to physical oceanographers such as Arakawa-C grids and Runge-Kutta time integration. We also avoid upwinding by using simple, centered finite differences for the advection terms.

- **Error control -** Embedded RK schemes are supported. Users can limit the velocity error using a variety of methods (local, PI, or PID controllers). From the user's perspective, all that is needed is an error tolerance specified in a text-based input file and SOMAR will take care of the rest.

- **Python post-processing -** Analyzing SOMAR's output can be accomplished in [VisIt](https://visit-dav.github.io/visit-website/index.html 'The VisIt website'), [ParaView](https://www.paraview.org/ 'The ParaView website'), or with simple Python scripts. Python can also be used to create initial conditions, create custom force functions, and process data on-the-fly.

- **Highly parallelizable -** SOMAR uses MPI and is build on the [Chombo framework.](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)


Future work
-----------
- **Immersed boundary method [in progress] -** Work is currently underway to incorporate immersed boundaries into SOMAR. This will allow us to incorporate complex boundaries provided by digital elevation models. This effort is largely complete and will become available once its associated manuscript is published.
- **Free surface [proposed] -** This would extend SOMAR's abilities to resolve fast-moving gravity waves at the upper boundary.


Documentation
-------------
Up-to-date documentation is available at [https://somar.readthedocs.io](https://somar.readthedocs.io).
