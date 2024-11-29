Quick start
===========

Preliminaries
-------------
To properly compile and run SOMAR, you will need

* A C++ >=17 compiler (tested with g++ and icpp).
* A fortran >= 77 compiler (tested with gfortran and ifort). Even if using ifort, gfortran must be present.
* Python >= 3.7m
* hdf5 library built with MPI support
* h5py linked to the hdf5 library with MPI support.
* mpi4py
* Perl >= 5.0
* lapack library
* An MPI environment (only if running in parallel)
* git

It is important to realize that the Intel compiler uses the std library headers provided by gcc.
This can be a problem if the gcc version the Intel compiler uses is too old to support C++17
features of the std library. Caveat emptor.

Get the software
----------------
SOMAR can be downloaded at its `GitHub repository <https://github.com/MUON-CFD/SOMAR>`_
or directly from the command line:

.. code-block:: console

    git clone https://github.com/MUON-CFD/SOMAR.git
    cd SOMAR
    git submodule init
    git submodule update





Configuring the Python environment and its dependencies
-------------------------------------------------------
We need a properly configured Python environment to compile and run SOMAR.
Especially on HPC
hardware, the provided default Python environment may have conflicts. It is thus recommended to run SOMAR in
a dedicated Python environment created with conda (or mamba). In the root directory of SOMAR there is a conda
environment file that can be used to rapidly create a conda environment that will work for SOMAR, using the GNU compiler.
It should work on a workstation and should get you going.
Conda is also highly recommended when working on a HPC environment, but in that case you probably want to use
the compiler/mpi enviroment provided by the HPC, to take advantage of the full capabilities of the hardware. You should inquire if
hdf5 libraries with MPI support are available, and use those. You should then use conda to create a barebone environment and use the
provided Python Script to install mpi4py and h5py.

For now, we proceed assuming that you are setting SOMAR on a workstation. Before proceeding, make sure that PYTHONPATH is unset. If using bash,

.. code-block:: console

    unset PYTHONPATH

and remove it from your .basrc if you happen to have it defined there.

Next, we create a conda environment using the provided environment file environment_lite.yml.
But first, edit the name block in environment_lite.yml if you plan to use a name other than SOMAR for the environment.
Also set the variable SOMAR_ROOT to point to the directory where you have downloaded SOMAR.

Now it is time to create the conda environment.

.. code-block:: console

    cd /path/to/SOMAR
    conda env create -f environment_lite.yml

This environment is based on Python 3.9. Note that it will install paraview, which you can
use to visualize outputs.


To activate the conda environment,

.. code-block:: console

    conda activate SOMAR

You will note that (SOMAR) is now displayed to the left of your prompt.

We need to copy the script Setup.sh so that it is called when the environment is initialized.

.. code-block:: console

    cp Setup.sh $CONDA_PREFIX/etc/conda/activate.d/
    cp Cleanup.sh $CONDA_PREFIX/etc/conda/deactivate.d/
    conda deactivate
    conda activate

Now you should have a working environment. You can test it as indicated below. If everything checks,
you are ready to compile SOMAR.

Working with existing hdf5 libraries and/or Intel compiler (HPC environment)
----------------------------------------------------------------------------
Here, we assume that your HPC environment is
managed with module. So first thing, make sure you have conda (or mamba), your compiler and your mpi module loaded.
In our example, we use gcc, mvapich and anaconda (to provide conda)

.. code-block:: console

    module list
    Currently Loaded Modules:
    1) gcc/9.1.0   2) mvapich2_2.3.3/gcc_9.1.0   3) anaconda/2019.10



Next, we create a basic conda environment

.. code-block:: console

conda create -n SOMAR_HPC -c conda-forge python=3.9
Note we have selected python 3.9 for our environment.
Also note that conda may ask you to execute a command to add a block
to your .basrc file (or equivalent if using a different shell). Follow
those instructions before proceeding further.

Now it is time to activate the Conda environment and purge the pip cache (just to be sure)

.. code-block:: console

    conda activate SOMAR_HPC
    pip cache purge



You will note that (SOMAR_HPC) is now displayed to the left of your prompt.
Enter the SOMAR root directory and run the provided Python script to install the required
Python modules and the hdf5 library (you may have a version of libhdf5 already on your system
but this ensures that hdf5 is compiled with the right mpi environment). By default, it uses mpicc.
For Intel, use mpiicc or for the latest version mpiicx.

.. code-block:: console


    cd SOMAR
    export MPICC=<your_mpicc>; PythonScripts/installRequiredModules.py

This will take some time. If everything goes well, the proper Python environment
should be set.

Testing the h5py/mpi4py environment
-----------------------------------

To make sure everything works

.. code-block:: console

    cd PythonScripts
    python testMPI.py

If everything is working, you should see


    My rank is 0

    My size is 1






Running the software without MPI
--------------------------------
To test drive the software on your machine, change into the :code:`execLockExchange`
folder and build the code in 2D, serial mode (remember to activate your conda environment first!)

.. code-block:: console

    cd exec/LockExchange
    ./buildall --D=2 --Serial

The resulting executable takes an input parameter, the input file. A sample
input file for a 2D run is provided. The output files will be placed in two
folders called :code:`hdf5_output` and :code:`check_points`. If these do not
exist in the :code:`exec/LockExchange` folder, create them before running the
simulation.

.. code-block:: console

    ./Somar_2D.Serial.gcc.ex inputs.basicTest2D.research1

If you wish to compile and run with the Intel compiler suite, use the
:code:`--IntelCompiler` switch. A

.. code-block:: console

    ./buildall --D=2 --Serial --IntelCompiler
    ./Somar_2D.Serial.intel.ex inputs.basicTest2D.research1


Running the software with MPI
-----------------------------
If you have MPI installed, then the compile-build process is similar.
For example, if your machine has 32 cores, you can compile using the
:code:`-j64` switch (2 jobs per core is usually fastest). Then, you can run
using the standard MPI commands. Note that the 2D lock exchange demo is set up
so that the domain neatly decomposes over 18 MPI ranks, so using all of the
available cores may not be ideal. Here, we compile with all available cores,
but only run on 18.

.. code-block:: console

    cd exec/LockExchange
    ./buildall --D=2
    mpirun -np 18 ./Somar_2D.MPI.gcc.ex inputs.basicTest2D.research1




Visualize the results
---------------------
The simulation will produce a number of .hdf5 files in the hdf5_output folder.
These can be opened with `VisIt <https://visit-dav.github.io/visit-website/>`_.
Within VisIt, click "Open" and select the entire series of :code:`*.hdf5` files.
Then, click `add -> pseudocolor -> b_total`. This is a density field initially
scaled to the interval :math:`[0,1]`.

If you double-click the b_total entry, the pseudocolor plot's properties dialog
will appear that will allow you to restrict pseudocolor's range. Choose the
minimum value to be 0 and the maximum to be 1.

Now, click draw. The initial lock exchange setup will appear with light water
(blue) on the left and heavy water (red) on the right. (If you don't like the
choice of colors, feel free to re-open the properies dialog and change them.)

There is a time slider that you can use to watch the flow evolve. By moving the
slider all the way to the right, you will render the output from the very last
timestep.

To view the grids, click `Add -> mesh`, then click `draw`. I suggest you open the
mesh plot's properties dialog and reduce the opacity to about 15\%. Now, you
can see the location of the fine grids. Alternatively, you can remove the mesh
plot and choose `Add -> subset -> levels`. Again, open the subset plot's
properties dialog and choose wireframe. You can also choose what colors
represent each AMR level. I suggest making level 0 (the coarsest level)
completely transparent and choosing a bright color like yellow or purple for
level 1 (the refined level). Click `draw`` and an outline of the fine grids will
appear that is a bit less invasive than the mesh plot.
