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
    subprocess.check_call([sys.executable, "-m", "pip",
                           "install", package,  "--upgrade"])
    

def runCommand(command):
    subprocess.check_call(command.split())

def main():
    import shutil
    """ installs what needs to be installed """
    HOME = os.path.expanduser('~')
    NeedSys = True
    NeedSite = True
    NeedInsert = True
    os.environ['CC'] = os.environ['MPICC'] if 'MPICC' in os.environ.keys() else 'mpicc'
    os.environ['HDF5_MPI'] = 'ON'
    os.environ['HDF5_DIR'] = os.getcwd()+'/lib/hdf5/'

# set the version of the packages needed. First we check if the user has defined
# a version, otherwise, we fall back on default
    MPI4PY_VER = "3.1.6"
    #H5PY_VER = os.environ['H5PY_VER'] if 'H5PY_VER' in os.environ.keys() else "3.2.1"
    HDF_VER = os.environ['HDF_VER'] if 'HDF_VER' in os.environ.keys() else "0"
    HDF_REL = os.environ['HDF_REL'] if 'HDF_REL' in os.environ.keys() else "1.12"
    #try:
    #    with open(HOME+'/.startup.py', 'r') as f:
    #        if 'import sys' in f.read():
    #            NeedSys = False
    #    with open(HOME+'/.startup.py', 'r') as f:
    #        if 'import site' in f.read():
    #            NeedSite = False
    #    with open(HOME+'/.startup.py', 'r') as f:
    #        if 'sys.path.insert(0,site.USER_SITE)' in f.read():
    #            NeedInsert = False
        # if we are here, .startup.py exists but does not have what we need
    #    with open(HOME+'/.startup.py', 'r') as f:
    #        original = f.read()
    #    prependLines = str()
    #    postpendLines = str()
    #    if NeedSys:
    #        prependLines += 'import sys \n'
    #    if NeedSite:
    #        prependLines += 'import site \n'
    #    if NeedInsert:
    #        postpendLines += 'sys.path.insert(0,site.USER_SITE) \n'
    #    with open(HOME+'/.startup.py', 'w') as f:
    #        f.write(prependLines+original+postpendLines)
    #except:
        # file does not exists
    #    prependLines = str()
    #    postpendLines = str()
    #    if NeedSys:
    #        prependLines += 'import sys \n'
    #    if NeedSite:
    #        prependLines += 'import site \n'
    #    if NeedInsert:
    #        postpendLines = 'sys.path.insert(0,site.USER_SITE) \n'
    #    with open(HOME+'/.startup.py', 'w') as f:
    #        f.write(prependLines+postpendLines)

    
    install('colorama')
    
    if shutil.which('gfortran') is None:
        from colorama import Fore, Back, Style 
        print(Style.BRIGHT + Fore.RED + "gfortran is needed by SCons, even if you plan to use Intel or CLang compilers")
        print("Install gfortran for you distro before proceeding!!!")
        print(Style.RESET_ALL)
        sys.exit(1)
    
    
    

    import tarfile
    install("wget")
    import wget
    
    import shutil

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

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-binary=mpi4py', 'mpi4py==3.1.6', '--ignore-installed', '--upgrade'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--no-binary=h5py', 'h5py ', '--ignore-installed',  '--upgrade'])
    
    #        runCommand('git clone https://github.com/h5py/h5py.git')
    #        os.chdir('h5py')
    #        install('.')
            #subprocess.check_call([sys.executable, '-m', 'pip', 'install', '.', '--user', '--upgrade'])
    #        os.chdir('../')
    #        runCommand('rm -rf h5py')


if __name__ == '__main__':
    main()
