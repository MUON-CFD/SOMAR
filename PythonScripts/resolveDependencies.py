#!/usr/bin/env python3
import subprocess
import pkg_resources
import os
import sys
import site
import shutil
# ensures that the local version of a package has precendence over the system packages
sys.path.insert(0, site.USER_SITE)


def install(package):
    print("Installing package "+package)
    subprocess.check_call([sys.executable, "-m", "pip",
                           "install", package, "--user", "--upgrade"])


def runCommand(command):
    subprocess.check_call(command.split())


if __name__ == '__main__':
    """ check if the dependencies are satisfied. If not, installs what needs to be installed """
    print(os.getcwd())
    HOME = os.path.expanduser('~')
    NeedSys = True
    NeedSite = True
    NeedInsert = True
    

# set the version of the packages needed. First we check if the user has defined
# a version, otherwise, we fall back on default
    MPI4PY_VER = os.environ['MPI4PY_VER'] if 'MPI4PY_VER' in os.environ.keys() else "3.0.3"
    H5PY_VER = os.environ['H5PY_VER'] if 'H5PY_VER' in os.environ.keys() else "3.2.1"
    HDF_VER = os.environ['HDF_VER'] if 'HDF_VER' in os.environ.keys() else "0"
    HDF_REL = os.environ['HDF_REL'] if 'HDF_REL' in os.environ.keys() else "1.12"
    try:
        with open(HOME+'/.startup.py', 'r') as f:
            if 'import sys' in f.read():
                NeedSys = False
        with open(HOME+'/.startup.py', 'r') as f:
            if 'import site' in f.read():
                NeedSite = False
        with open(HOME+'/.startup.py', 'r') as f:
            if 'sys.path.insert(0,site.USER_SITE)' in f.read():
                NeedInsert = False
        # if we are here, .startup.py exists but does not have what we need
        with open(HOME+'/.startup.py', 'r') as f:
            original = f.read()
        prependLines = str()
        postpendLines = str()
        if NeedSys:
            prependLines += 'import sys \n'
        if NeedSite:
            prependLines += 'import site \n'
        if NeedInsert:
            postpendLines += 'sys.path.insert(0,site.USER_SITE) \n'
        with open(HOME+'/.startup.py', 'w') as f:
            f.write(prependLines+original+postpendLines)
    except:
        # file does not exists
        prependLines = str()
        postpendLines = str()
        if NeedSys:
            prependLines += 'import sys \n'
        if NeedSite:
            prependLines += 'import site \n'
        if NeedInsert:
            postpendLines = 'sys.path.insert(0,site.USER_SITE) \n'
        with open(HOME+'/.startup.py', 'w') as f:
            f.write(prependLines+postpendLines)

    try:
        pkg_resources.require('colorama')
    except:
        install('colorama')
    
    if shutil.which('gfortran') is None:
        from colorama import Fore, Back, Style 
        print(Style.BRIGHT + Fore.RED + "gfortran is needed by SCons, even if you plan to use Intel or CLang compilers")
        print("Install gfortran for you distro before proceeding!!!")
        print(Style.RESET_ALL)
        sys.exit(1)
    
    try:
        pkg_resources.require('mpi4py>='+MPI4PY_VER)
    except:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-binary=mpi4py',
                               'mpi4py=='+MPI4PY_VER, '--ignore-installed', '--user', '--upgrade'])

    try:
        pkg_resources.require('h5py>='+H5PY_VER)
        # print("checking form h5py>="+H5PY_VER)
        # import h5py as h5
        # LOCAL_VER=h5.__version__.split('.')
        # for s,p in zip(LOCAL_VER,H5PY_VER.split('.')):
        #     if int(s) < int(p):
        #         raise ValueError(' ')

        from mpi4py import MPI
        import h5py as h5
        # check for parallel
        FileHandle = h5.File(
            'test.h5', 'a' if os.path.exists('test.h5') else 'w', driver='mpio', comm=MPI.COMM_WORLD)
        FileHandle.close()
        runCommand('rm test.h5')
    except:
        print("h5py not found.  I will try to install it now. This may take a while...")
        print("as in 'time to grab a coffee' while")
        try:
            import tarfile
        except:
            install("tarfile")
            import tarfile
        try:
            import wget
        except:
            install("wget")
            import wget
        try:
            pkg_resources.require('cython>=0.29')
        except:
            install("cython")

        try:  # before recompiling, let's make sure we have what we need already at hand
            if not os.path.isfile('lib/hdf5/lib/libhdf5.settings'):
                print('hdf5 not found')
                raise ValueError(' ')
            with open('lib/hdf5/lib/libhdf5.settings') as f:
                if not 'HDF5 Version: '+HDF_REL+'.'+HDF_VER in f.read():
                    print('wrong version of hdf5 installed')
                    raise ValueError(' ')
            with open('lib/hdf5/lib/libhdf5.settings') as f:
                if not 'Parallel HDF5: yes' in f.read():
                    print('hdf5 was not compiled with parallel support')
                    raise ValueError(' ')
            with open('lib/hdf5/lib/libhdf5.settings') as f:
                if not 'Libraries: static, shared' in f.read():
                    print('hdf5 was not compiled with shared libraries support')
                    raise ValueError(' ')
        except:
            import shutil
            os.environ['CC'] = 'mpicc'
            HDF5 = 'hdf5-'+HDF_REL+'.'+HDF_VER
            url = 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-' + \
                HDF_REL+'/'+HDF5+'/src/'+HDF5+'.tar.gz'
            filename = wget.download(url)
            tarfile.open(filename).extractall()
            runCommand('rm '+HDF5+'.tar.gz')
            os.chdir(HDF5)
            runCommand('./configure --enable-silent-rules --enable-tests=no --enable-tools=no --enable-parallel --enable-shared --prefix='+os.getcwd()+'/../lib/hdf5')
            runCommand('make -j 16')
            runCommand('make install')
            runCommand('make clean')
            os.chdir('../')
            runCommand('rm -r '+HDF5)

        os.environ['HDF5_MPI'] = 'ON'
        os.environ['HDF5_DIR'] = os.getcwd()+'/lib/hdf5/'
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-binary=h5py', 'h5py=='+H5PY_VER, '--ignore-installed', '--user', '--upgrade'])
            raise ValueError(' Installing from distribution failed, trying compiling from sources')
        except:
            runCommand('git clone https://github.com/h5py/h5py.git')
            os.chdir('h5py')
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', '.', '--user', '--upgrade'])
            os.chdir('../')
            runCommand('rm -rf h5py')
