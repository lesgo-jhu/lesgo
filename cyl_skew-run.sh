#!/bin/bash
NPROC=$1
EXEC=cylinder_skew
TAIL_NUM=100
USE_GDB=YES
if [ -n "$NPROC" ]; then
  /bin/true
else
  echo Error: please specify number of processors.
  exit 1
fi
if [ $USE_GDB == "YES" ]; then
  mpiexec -gdb -np $NPROC ./$EXEC 
else
  mpiexec -np $NPROC ./$EXEC >& $EXEC.out < /dev/null &
  tail -n $TAIL_NUM -f $EXEC.out  
fi
exit 0

