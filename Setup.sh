#/bin/sh 

python_version=$(python --version 2>&1)
version=$(echo $python_version | awk '{print $2}')
major_minor=$(echo $version | cut -d '.' -f 1,2)

export PYTHON_VER=$major_minor
export PYTHON_LIB_PATH=$CONDA_PREFIX/lib/
export PYTHON_INCLUDE_PATH=$CONDA_PREFIX/include/python$major_minor
export LD_LIBRARY_PATH_OLD=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$PYTHON_LIB_PATH:$LD_LIBRARY_PATH

