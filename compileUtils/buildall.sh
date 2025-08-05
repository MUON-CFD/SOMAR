#!/bin/bash
ARGS=${@:1}
EXEC_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd ../..
python3 ./SCons/scons.py $ARGS
cd $EXEC_PATH
python3 ../../SCons/scons.py $ARGS

