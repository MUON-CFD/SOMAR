#*******************************************************************************
#  SOMAR - Stratified Ocean Model with Adaptive Refinement
#  Developed by Ed Santilli & Alberto Scotti
#  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
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
#  https://github.com/MUON-CFD/SOMAR.
# ******************************************************************************/

import os
import SCons as SCons
import sys
import glob
import subprocess

print("using SCons version ", SCons.__version__)
Flags=processFlags()
# first build checks for dependencies
try: # only run if called from root directory.
     #subprocess.check_call([sys.executable, 'PythonScripts/resolveDependencies.py'])
     pass
except:
     pass



env=Environment(**setEnvOptions(Flags), ENV=os.environ, tools=['default', ADD_TOOLS] )

libDir=buildName('lib/', Flags)
# since the Include files are referenced to the src dir,
# the _F.H header are referenced there and will be build there
# this is how it was with GNUMakefile anyway. It is actually
# not a bad idea, it helps the IDE.
Progress(['-\r', '\\\r', '|\r', '/\r'], interval=5)
ChFHeader=[]
ChFNodes=variantglob(env, '*.ChF', recursive=True)
for file in ChFNodes:
     ChFHeader+=env._H(source=file)


# now we switch to the build dir
buildDir=buildName('build/', Flags)

env.Replace(VARIANT_DIR= buildDir)
VariantDir(variant_dir=env['VARIANT_DIR'],
           src_dir=env['SOURCE_DIR'], duplicate=0)


# env.Object can work on lists


ChFNodes=variantglob(env, '*.ChF', recursive=True)
CppNodes=variantglob(env, '*.cpp', recursive=True)


if Flags.StaticLib:
    ChFobj=env.Object(ChFNodes)
    Cppobj=env.Object(CppNodes)
    env.NoClean(env.Library(target=libDir+'/SOMAR', source=Cppobj+ChFobj))
else:
    ChFobj=env.SharedObject(ChFNodes)
    Cppobj=env.SharedObject(CppNodes)
    env.NoClean(env.SharedLibrary(target=libDir+'/SOMAR', source=Cppobj+ChFobj, SHLIBVERSION='1.0.99'))








