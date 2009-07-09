#!/bin/bash
NPROC=$1
EXEC=lesgo
TAIL_NUM=100
if [ -n "$NPROC" ]; then
  /bin/true
else
  echo Error: please specify number of processors.
  exit 1
fi
mpiexec -np $NPROC ./$EXEC >& $EXEC.out < /dev/null &
tail -n $TAIL_NUM -f $EXEC.out 
exit 0

