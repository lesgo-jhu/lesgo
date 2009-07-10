#!/bin/bash
NPROC=$1
EXEC=lesgo
TAIL_NUM=100
USE_GDB=no
if [ -n "$NPROC" ]; then
  /bin/true
else
  echo Error: please specify number of processors.
  exit 1
fi
if [ $USE_GDB == "yes" ]; then
  mpiexec -gdb -np $NPROC ./$EXEC 
else
  mpiexec -np $NPROC ./$EXEC >& $EXEC.out < /dev/null &
  tail -n $TAIL_NUM -f $EXEC.out  
fi
exit 0

