#! /usr/bin/env python3
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
#  https://github.com/MUON-CFD/somar.
# ******************************************************************************/

# this collection of function is used by SCons to build SOMAR
import os
import sys
import site

if site.USER_SITE in sys.path:
    sys.path.remove(site.USER_SITE)
sys.path.insert(0,site.USER_SITE) # ensures that the local version of a package has precendence over the system packages

#import pkg_resources
import SCons as SCons
import glob
import subprocess
import shutil


def ADD_TOOLS(env):
    """ Fortran and C++ files with extensions .f and .cpp are handled by SCons
    automatically. The conversion of Chombo Fortran files to .f files and _F.H files
    requires that we specify the recipe, which is done here.  The process of going from .ChF to .f file
    is done in three steps. First, the .ChF is processed by a perl script (actually,
    a collection) that substitutes the CH macros. The result is written to a file
    with extension .C. This is in turn passed to the C preprocessor which eliminates
    the #includes and sets further macros (e.g., REAL_T gets substituted by real*8
     real*4 depending on the compiler flag), and finally it is passed to a perl script that formats the output in fixed format 72 columns.
    If a user ever needs to add more recipes (e.g., on how to handle .cu files)
    should follow a similar template. For more info, see Builders in the SCons
    documentation.
    """

    compileUtils=env['compileUtils'] # the directory containing the Perl cripts
    IncludePath=env['IncludePath']
    DSwitches=env['DSwitches']
    Dim=env['Dim']

    action_ChfToC="perl -I "+compileUtils+"chfpp "+compileUtils+str("chfpp/uber.pl -f $SOURCE -p /dev/stdout -c /dev/null -D")+str(Dim)
    action_f= action_ChfToC+' | '+"g++ -P -E " + DSwitches  + " " + IncludePath + " - | perl -I "+compileUtils+str("chf\
pp ")+compileUtils+str("chfpp/fort72 > $TARGET")

    bld_Chf=Builder(suffix='.f', src_suffix='.ChF', action=action_f)

    action_ChfToH="perl -I "+compileUtils+"chfpp "+compileUtils+str("chfpp/uber.pl -f $SOURCE -p /dev/null -c $TARGET  -D")+str(Dim)
    bld_H=Builder(suffix='_F.H', src_suffix='.ChF', action=action_ChfToH)





    env.Append(BUILDERS={'CHF' : bld_Chf})
    env.Append(BUILDERS={'_H' : bld_H})
    from SCons import SCons
    static_obj, shared_obj = SCons.Tool.createObjBuilders(env)
    static_obj.add_src_builder('CHF')

    shared_obj.add_src_builder('CHF')



def variantglob(env, pattern, ondisk=True, source=True, strings=False,
                recursive=False):
    """ build a list of files that need to be processed traversing the 'SOURCE_DIR' tree
       but building relative to the 'VARIANT_DIR' tree.  """
    matches = []
    for root, dirs, filenames in os.walk(env['SOURCE_DIR'], followlinks=True):
        cwd = Dir(os.path.join(env['VARIANT_DIR'],
                               os.path.relpath(root, env['SOURCE_DIR'])))
        matches.extend(cwd.glob(pattern, ondisk, source, strings))
    return matches



class processFlags:
    """ Reads command line options and stores them.
    If the user needs to add flags, follow this template and
    add it to self """
    def __init__(self):
        AddOption('--D',
                  dest='Dim',
                  type='int',
                  nargs=1,
                  action='store',
                  metavar='DIMENSION',
                  default=2,
                  help='Dimension of build: must be 2 or 3')
        AddOption('--Debug',
                  dest='Debug',
                  action='store_true',
                  metavar='DEBUG',
                  help='build with debug options')
        AddOption('--StaticLib',
                    dest='StaticLib',
                    action='store_true',
                    metavar='STATIC_LIB',
                    help='build a static version of the library')
        AddOption('--Serial',
                dest='MPI',
                action='store_false',
                metavar='SERIAL',
                help='build the serial version of the code')
        AddOption('--Profile',
                dest='Profile',
                action='store_true',
                metavar='PROFILE',
                help='turn on the -pg option for code profiling')
        AddOption('--OpenMP',
                    dest='OpenMP',
                    action='store_true',
                    metavar='OpenMP',
                    help='turn on the options to compile with OpenMP support')
        AddOption('--noPython',
                    dest='noPython',
                    action='store_true',
                    metavar='Python',
                    help='Will compile without support for python interpreter. Warning: no plot nor checkpoint files will be created')
        AddOption('--IntelCompiler',
                    dest='IntelC',
                    action='store_true',
                    metavar='Intel',
                    help='Will use the Intel compiler. If either ifort or ipcp are not present, it will revert to using gcc')
        AddOption('--Clang',
                    dest='Clang',
                    action='store_true',
                    metavar='clang',
                    help='Will use the clang compiler for c++ stuff. If clang++ is not present, it will revert to using gcc. Fortran files are still handled by gfortran')
        AddOption('--Fastidious',
                    dest='Fastidious',
                    action='store_true',
                    metavar='fastidious',
                    help='Turns on extra compiler warnings')

        self.Dim=GetOption('Dim')
        if (self.Dim>3) or (self.Dim<2):
            raise ValueError("Dim must be either 2 or 3")
        self.MPI=GetOption('MPI') if GetOption('MPI') is not None else True
        self.Debug=GetOption('Debug') if GetOption('Debug') is not None else False
        self.StaticLib=GetOption('StaticLib') if GetOption('StaticLib') is not None else False
        self.Profile=GetOption('Profile') if GetOption('Profile') is not None else False
        self.OpenMP=GetOption('OpenMP') if GetOption('OpenMP') is not None else False
        self.noPython=GetOption('noPython') if GetOption('noPython') is not None else False
        self.IntelCompiler=GetOption('IntelC') if GetOption('IntelC') is not None else False
        self.ClangCompiler=GetOption('Clang') if GetOption('Clang') is not None else False
        self.Fastidious = GetOption('Fastidious') if GetOption('Fastidious') is not None else False

        # check for the intel compiler to be actually present
        if self.IntelCompiler:
            self.IntelCompiler= (shutil.which('ifort') is not None) and (shutil.which('icx') is not None)
            if not self.IntelCompiler: print("Intel compiler missing, reverting to gcc")
        # check of the clang compiler to be actually present
        if self.ClangCompiler:
            self.ClangCompiler=(shutil.which('clang++') is not None)
            if not self.ClangCompiler: print("clang compiler missing, reverting to gcc")

        if self.ClangCompiler and self.IntelCompiler:
            raise ValueError("You must specify either --Clang or --IntelCompiler, not both.")

        self.GNUCompiler= not (self.ClangCompiler or self.IntelCompiler)

def setEnvOptions(Flags, src_dir='src', root_dir='.'):
    """ Used to set compiler flags and switches. CXXFLAGS (a list) are passed
    to the C++ compiler when it generates .o files. CPPSwitches are passed to the
    C preprocessor prepending -D in front of each. Note that CPPSwitches is a dictionary.
    If a switch does not need a value (e.g., because it is used as #ifdef MySwitch) the value should be set to None.
     If a switch needs a value (e.g., Dim) it should be set as the value for that particular key.
     F77FLAGS (a list) are passed to
    the Fortran compiler which works on the .f files created from Chombo Fortran files. """
    import os
    # flags and switches that are common to all builds
    CXXFLAGS=['-std=c++17'
              ,'-Wall'                  # reasonable and standard
			#   ,'-Wextra'                # reasonable and standard
			#   ,'-Wshadow'               # warn the user if a variable declaration shadows one from a parent context
			#   ,'-Wnon-virtual-dtor'     # warn the user if a class with virtual functions has a non-virtual destructor. This helps catch hard to track down memory errors
			#   ,'-Wduplicated-cond'      # (only in GCC >= 6.0) warn if if / else chain has duplicated conditions
			#   ,'-Wduplicated-branches'  # (only in GCC >= 7.0) warn if if / else branches have duplicated code
			#   ,'-Wlogical-op'           # (only in GCC) warn about logical operations being used where bitwise were probably wanted
			#   ,'-Wnull-dereference'     # (only in GCC >= 6.0) warn if a null dereference is detected
			# #   ,'-Wlifetime'             # (only special branch of Clang currently) shows object lifetime issues
			#   ,'-Wpedantic'             # (all versions of GCC, Clang >= 3.2) warn if non-standard C++ is used
			#   ,'-pedantic'              # Warn on language extensions
			#   ,'-Wcast-align'           # warn for potential performance problem casts
			#   ,'-Wunused'               # warn on anything being unused
			#   ,'-Woverloaded-virtual'   # warn if you overload (not override) a virtual function
			#   ,'-Wconversion'           # warn on type conversions that may lose data
			#   ,'-Wsign-conversion'      # (Clang all versions, GCC >= 4.3) warn on sign conversions
			#   ,'-Wuseless-cast'         # (only in GCC >= 4.8) warn if you perform a cast to the same type
			#   ,'-Wdouble-promotion'     # (GCC >= 4.6, Clang >= 3.8) warn if float is implicit promoted to double
			#   ,'-Wformat=2'             # warn on security issues around functions that format output (ie printf)
            #   ,'-fpermissive'
            # , '-m64'
			#   ,'-DCH_NTIMER'
              ]
    CPPSwitches={'CH_SPACEDIM' : Flags.Dim,
            'CH_Linux' : None,
            'CH_USE_COMPLEX' : None,
            'CH_USE_64' : None,
            'CH_USE_DOUBLE' : None,
            'CH_USE_PYTHON' : None,
            'CH_FORT_UNDERSCORE' : None,
            'CH_LANG_CC' : None,
			'CH_NTIMER' : None
            }
    F77FLAGS=[#'-m64'
             ]
    # Here we customize to specific builds. Note that we append (+=) to list and
    # define new entries for the dictionary of switches.
    if Flags.Fastidious:
        CXXFLAGS+=['-Wextra'                # reasonable and standard
		          ,'-Wshadow'               # warn the user if a variable declaration shadows one from a parent context
			      ,'-Wnon-virtual-dtor'     # warn the user if a class with virtual functions has a non-virtual destructor. This helps catch hard to track down memory errors
        ]
        if not ( Flags.IntelCompiler or Flags.ClangCompiler):
            CXXFLAGS+=['-Wduplicated-cond'      # (only in GCC >= 6.0) warn if if / else chain has duplicated conditions
			          ,'-Wduplicated-branches'  # (only in GCC >= 7.0) warn if if / else branches have duplicated code
			          ,'-Wlogical-op'           # (only in GCC) warn about logical operations being used where bitwise were probably wanted
			          ,'-Wnull-dereference'     # (only in GCC >= 6.0) warn if a null dereference is detected
                      ,'-Wno-unused-parameter'  # This warning is just annoying.
			]
    if Flags.Debug:
        if Flags.IntelCompiler or Flags.ClangCompiler:
            CXXFLAGS+=['-g', '-O0']
            F77FLAGS+=['-g', '-O0']
        else:
            CXXFLAGS+=['-ggdb3', '-O0']
            F77FLAGS+=['-ggdb3', '-O0', '-fbounds-check']
        CPPSwitches['DEBUG']=None
        CPPSwitches['CH_USE_SETVAL']=None # this is a bit unfortunate, but it means that it sets corresponding switch
    else:
        CXXFLAGS+=['-O3']
        F77FLAGS+=['-O3', '-funroll-loops']
        CPPSwitches['NDEBUG']=None # see comment above

    if Flags.noPython: del CPPSwitches['CH_USE_PYTHON']
    if Flags.Profile:
        CXXFLAGS+=['-pg', '-g']
        F77FLAGS+=['-pg', '-g']

        # # For Linaro MAP
        # CXXFLAGS+=['-g1', '-O3', '-fno-inline', '-fno-optimize-sibling-calls']
        # F77FLAGS+=['-g1', '-O3', '-fno-inline', '-fno-optimize-sibling-calls']

        # CXXFLAGS+=['-g', '-O2', '-shared-intel', '-debug inline-debug-info', '-D TBB_USE_THREADING_TOOLS', '-parallel-source-info=2']
        # F77FLAGS+=['-g', '-O2', '-funroll-loops', '--info-for-profiling']

    if Flags.MPI:
        CPPSwitches['CH_MPI']=None

    if Flags.OpenMP:
        if Flags.IntelCompiler or Flags.ClangCompiler:
            raise ValueError("OpenMP not supported yet for intel or clang")
        F77FLAGS+=['-fopenmp']
        CXXFLAGS+=['-fopenmp']

    if Flags.IntelCompiler:
        CXXFLAGS+=['-diag-disable=10441']    # The Intel(R) C++ Compiler Classic (ICC) is deprecated and will be removed from product release in the second half of 2023. The Intel(R) oneAPI DPC++/C++ Compiler (ICX) is the recommended compiler moving forward. Please transition to use this compiler.

    # finally we go about the business of creating a list of subdirectories to inspect for compile files.
    compileUtils=root_dir+'/compileUtils/'
    highLevelDirs=['Chombo', 'Basics', 'AnisotropicChombo', 'Calculus', 'IB', 'SOMAR']
    IncDirs=[] # will contain all the subdirs in the source tree
    srcRoot=root_dir+'/src/'
    if Flags.noPython:
        excludeDirs = set([os.getcwd()+ '/'+ srcRoot + 'Grade0_Chombo/PyGlue'])
    else:
        excludeDirs = set([])
    for n,highLevelDir in enumerate(highLevelDirs):
        IncDirs.append([x[0] for x in os.walk(os.getcwd()+ '/'+ srcRoot + 'Grade' + str(n) + '_' + highLevelDir + '/',followlinks=True) if x[0] not in excludeDirs])
    IncludePath=str() # this is passed to the perl scripts.
    for IncDir in IncDirs:
        for dir in IncDir:
            IncludePath+= str("-I")+dir+str(" ")


# the following is a bit of a hack but until I figure out how to tell scons that the .C files are to preprocessed
# by cpp it will have to do.
    DSwitches=str()
    for key,value in CPPSwitches.items():
        if key == 'CH_SPACEDIM':
            DSwitches+=str("-D")+key+str("=")+str(value)+str(" ")
        elif key == 'CH_LANG_CC':
            DSwitches+='-DCH_LANG_FORT '
        else:
            DSwitches+=str("-D")+key+str(" ")

    # the location of libpythonX.Ym and Python.h are obtained from environmental variables
    # set by the user. If not set, we will try to look in the usual places, but do not expect a
    # miracle... If compilation files complaining that Python.H cannot be found, it means
    # that either PYTHON_INCLUDE_PATH was not set and the default fails, or it was set incorrectly
    # ditto for ld cannot find libpython3.6m. On Ubuntu, PythonLibPath does not seem to be needed,
    # the linker knows where to find it, but I have added here anyway for clarity.
    # The version of Python to use can be specified setting the environmental variable PYTHON_VER
    # If not specified, the default value 3.7m is used.

    try:
        PYTHON_VER=os.environ['PYTHON_VER']
    except:
        PYTHON_VER="3.7m"
    if Flags.noPython:
        PythonLibPath=None
        PythonIncludePath=None
    else:
        try:
            PythonLibPath=os.environ['PYTHON_LIB_PATH']
        except:
            PythonLibPath = '/usr/lib/x86_64-linux-gnu' if os.path.isdir('/usr/lib/x86_64-linux-gnu') else None

        try:
            PythonIncludePath=os.environ['PYTHON_INCLUDE_PATH']
        except:
            PythonIncludePath= '/usr/include/python'+PYTHON_VER if os.path.isdir('/usr/include/python'+PYTHON_VER) else None

    # finally all the information is packaged into a dictionary which will be passed to SCons
    env_options = {
        "CPPPATH" : IncDirs + [PythonIncludePath or str('')],
        "LIBPATH" : ['.',root_dir+buildName('/lib/',Flags), PythonLibPath or str('')],
        "FORTRANFLAGS" : F77FLAGS,
        "SHF77FLAGS" : F77FLAGS,
        "CPPDEFINES" : CPPSwitches,
        "CXXFLAGS" : CXXFLAGS,
        # "LINKFLAGS" : ['-pg -diag-disable=10441'] if Flags.Profile else ['-diag-disable=10441'],
        "SOURCE_DIR" : src_dir,
        "VARIANT_DIR" : src_dir
    }

    # here we sepcify the name of the compiler.
    if Flags.MPI:
        env_options['CXX']='mpiicpc' if shutil.which('mpiicpc') is not None and Flags.IntelCompiler else 'mpic++'

    else:
        if Flags.IntelCompiler:
            env_options['CXX']= 'icx'
        elif Flags.ClangCompiler:
            env_options['CXX']= 'clang++'
        else:
            env_options['CXX']= 'g++'


    env_options['FORTRAN']='ifort' if Flags.IntelCompiler else 'gfortran'

    #env_options['CXX'] += ' -pg' if Flags.Profile else []
    env_options['IncludePath']=IncludePath
    env_options['DSwitches']=DSwitches
    env_options['compileUtils']=compileUtils
    env_options['Dim']=Flags.Dim

    return env_options

def buildName(prefix,Flags):
    """ append to prefix the string with dim, MPI or Serial and Debug. """
    prefix+=str(Flags.Dim)+'D'
    prefix+= '.MPI' if Flags.MPI else '.Serial'
    if Flags.Debug: prefix+='.debug'
    if Flags.Profile: prefix+='.profile'
    if Flags.OpenMP: prefix+='.OpenMP'
    if Flags.noPython: prefix+='.noPython'
    if Flags.StaticLib: prefix+='.static'
    if Flags.IntelCompiler: prefix+='.intel'
    if Flags.ClangCompiler: prefix+='.clang'
    if Flags.GNUCompiler: prefix+='.gcc'

    return prefix
