#*******************************************************************************
#  SOMAR - Stratified Ocean Model with Adaptive Refinement
#  Developed by Ed Santilli & Alberto Scotti
#  Copyright (C) 2020
#    Jefferson (Philadelphia University + Thomas Jefferson University) and
#    University of North Carolina at Chapel Hill
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
#  USA
#
#  For up-to-date contact information, please visit the repository homepage,
#  https://github.com/somarhub.
# ******************************************************************************/

import os
import SCons as SCons
import sys
import glob
import colorama
from colorama import Fore, Style
print("using SCons version ", SCons.__version__)

try:
    PYTHON_VER = os.environ['PYTHON_VER']
except:
    PYTHON_VER = "3.7m"

Flags=processFlags()



env=Environment(**setEnvOptions(Flags, src_dir='.', root_dir='../..'), ENV=os.environ, tools=['default', ADD_TOOLS] )
rootLibDir=os.getcwd()+'/../../lib/'
SomarLibDir=buildName(rootLibDir, Flags)

env['RPATH']=SomarLibDir



# We are ready to compile!

# set the name of the executable and the libPath
execName=buildName('Somar_',Flags)
execName+= '.ex'

# uncomment the following lines and add paths to specific libraries
#libPath=
#env['LIBPATH'].append(libPath)




# in case we forgot about these
if not os.path.exists('check_points'):
    os.makedirs('check_points')
if not os.path.exists('hdf5_output'):
    os.makedirs('hdf5_output')



# create an internal list of header files that need to be generated
# to deal with Chombo Fortran files.
ChFNodes=variantglob(env, '*.ChF', recursive=True)
ChFHeader=[]
for file in ChFNodes:
    ChFHeader.append(env._H(file))

# now we switch to the build dir to assemble the list of files
# to be converted to object files
currentDirName=os.getcwd().split(os.sep)[-1]
buildDir=buildName('../../build/'+currentDirName+'.',Flags)
env.Replace(VARIANT_DIR= buildDir)
VariantDir(variant_dir=env['VARIANT_DIR'],
           src_dir=env['SOURCE_DIR'], duplicate=0)
ChFNodes=variantglob(env, '*.ChF', recursive=True)
CppNodes=variantglob(env, '*.cpp', recursive=True)
# the libraries that need to be linked with.
SomarLib=['SOMAR']
# add here if you need other libraries, but careful, the order may matter
OtherLibs=[]
if not Flags.noPython:
    OtherLibs.append('python'+PYTHON_VER)

OtherLibs.append('lapack')
OtherLibs.append('m')
if Flags.OpenMP:
    OtherLibs.append('gomp')

if Flags.StaticLib or Flags.Debug:
    OtherLibs.append('gfortran')

if Flags.IntelCompiler:
    OtherLibs.append('ifcore')
    OtherLibs.append('ifcoremt')
    OtherLibs.append('pthread')

# make the object files
ChFObj=env.Object(ChFNodes) # Chombo Fortran
CppObj=env.Object(CppNodes) # C++ files
SomarObj=env.Object(source='../somar.cpp', target='./somar.o') # main()

# link into the executable
env.Program(execName, source=SomarObj+ChFObj+CppObj, LIBS=SomarLib+OtherLibs)

print(" executable name is " + Style.BRIGHT+Fore.RED+'\033[1m '+execName+'\033[0m')
print(Style.RESET_ALL)

