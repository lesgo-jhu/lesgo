#!/usr/bin/bash
#@output=out.$(jobid)
#@error=err
#@job_type=parallel
#@node=1
#@total_tasks=8
#@node_usage=shared
#@class=share
#@queue

export MP_SHARED_MEMORY=yes
export OBJECT_MODE=64
export XLSMPOPTS="parthds=8:spins=500000:yields=50000"

date
./interp
date

