#!/usr/bin/bash
suffix=$(date +%T--%m-%d-%Y)
cd output
list=$(ls -x *-avg_stats.dat)
echo "suffix:==${suffix}=="
echo "list:==${list}=="
for file in $list; do
  echo "file:==${file}=="
  cp $file $file.$suffix
done
