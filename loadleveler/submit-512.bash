#!/usr/bin/bash
#@output=out.$(jobid)
#@error=err
#@job_type=parallel
#@network.MPI=csss,shared,us
#@node=2
#@total_tasks=64
#@node_usage=not_shared
#@class=com_rg32
#@queue

export MP_SHARED_MEMORY=yes

#--note after a series of runs, only the last output will be left
#echo "removing velXXXXXX.out.c* files"
#rm -f ./output/vel[0-9][0-9][0-9][0-9][0-9][0-9].out.c*

date
./les.go
date

RESUBMIT='FALSE'

if [ -e ./resubmit ]; then
  N=$(cat ./resubmit)
  echo "file resubmit exists and requests $N more job submissions"
  if [[ $N > 0 ]]; then
    RESUBMIT='TRUE'
    (( N-- ))
    echo $N > ./resubmit
  fi
fi

if [ $RESUBMIT == 'TRUE' ]; then
  echo "$0: resubmitting job"
  date
  sleep 120
  date
  /ssg/loadl/bin/llsubmit $0
else
  echo "$0: not resubmitting job"
fi
