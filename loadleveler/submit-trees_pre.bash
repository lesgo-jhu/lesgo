#!/usr/bin/bash
#@output=out.$(jobid)
#@error=err
#@job_type=parallel
#@node=1
#@tasks_per_node=32
#@node_usage=not_shared
#@class=com_rg32
#@queue

ntot=64  #--this is npus to main code
ncpu=32  #--this should be ncpus for this script (not real simulation)
m=$((ntot/ncpu))

for i in $(seq 0 $(($ncpu-1))); do
  cat > temp$i.dat << EOF
$(($m*$i))
$(($m*($i+1)-1))
EOF
done

date
for i in $(seq 0 $(($ncpu-1))); do
  ./trees_pre_ls < temp$i.dat &> out.$i &
done
wait
date
