#!/usr/bin/bash
#@output=out.$(jobid)
#@error=err
#@job_type=parallel
#@network.MPI=csss,shared,us
#@node=4
#@total_tasks=8
#@node_usage=shared
#@class=share
#@queue

export MP_SHARED_MEMORY=yes

date
./les.go
date
./save-avg_stats.bash

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
