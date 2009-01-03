#!/usr/bin/ksh
#--LSF commands
#BSUB -P 35151020
#BSUB -n 1
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -J serial
#BSUB -q share
#BSUB -W 3:00

date
./les.go
date

RESUBMIT='FALSE'

if [ -e ./resubmit ]; then
  N=$(cat ./resubmit)
  echo "file resubmit exists and requests $N more job submissions"
  if [[ $N > 0 ]]; then
    RESUBMIT='TRUE'
    N=$(( N-1 ))
    echo $N > ./resubmit
  fi
fi

if [ $RESUBMIT == 'TRUE' ]; then
  echo "$0: resubmitting job"
  date
  sleep 10
  date
  bsub < $0
else
  echo "$0: not resubmitting job"
fi
