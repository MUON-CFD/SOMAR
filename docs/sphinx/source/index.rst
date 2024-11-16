.. SOMAR documentation master file, created by
   sphinx-quickstart on Wed Jul 14 12:02:27 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Stratified Ocean Model with Adaptive Refinement (SOMAR)
===========================================================

.. note::
   This project is under active development.


SOMAR is a Boussinesq ocean model with several key features:

* Simulates 2D or 3D, non-hydrostatic flows.
* Can adaptively refine in space and time.
* Comes with a Large Eddy Simulation (LES).
* Can handle any grid stretching functions of the form :math:`\left( x, y, z \right) = \left( f(\xi), g(\eta), h(\zeta) \right)`.
* Can include a static, vertical background stratifications of temperature and salinity.
* Uses the familiar Arakawa-C grids and Runge-Kutta time integration.
* Highly parallelizable using `MPI <https://www.mpich.org/>`_ and the `Chombo framework <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_.
* Supports immersed boundaries (IB)


SOMAR is `free software <https://www.gnu.org/licenses/lgpl-2.1.html>`_ provided
jointly by `Thomas Jefferson University's College of Humanities and Sciences <https://www.jefferson.edu/academics/colleges-schools-institutes/humanities-sciences.html>`_
and `Arizona State University's School of Engineering <https://semte.engineering.asu.edu/>`_.


.. toctree::
    :caption: Tutorials
    :maxdepth: 2
    :hidden:

    quick_start
    eom
    new_sim
    ib
    amr_primer
    new_amr_sim
    refs

.. .. toctree::
..     :caption: Example simulations
..     :maxdepth: 2
..     :hidden:

..     Lock exchange
..     DJL solitary wave
..     ???

.. toctree::
    :caption: Software reference
    :maxdepth: 2
    :hidden:

    grade0_ref
    grade1_ref
    grade2_ref
    grade3_ref
    grade4_ref
    grade5_ref



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
