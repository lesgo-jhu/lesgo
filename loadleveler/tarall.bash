#!/bin/bash

#--vel part
nfiles=64
fprefix="vel.out.c"

for i in $(seq $((($nfiles+10)/10))); do

  start=$((10 * ($i - 1)))
  stop=$(($start + 9))
  sufs=$(seq $start $stop)

  files=""
  for s in $sufs; do
    if [ -e "$fprefix$s" ]; then
      files="$files $fprefix$s"
    fi
  done

  tar -cvf - $files | gzip > vel$i.tgz
  rm -f $files

done

#--phi part
tar -cvf - phi.out.c* | gzip > phi.tgz

