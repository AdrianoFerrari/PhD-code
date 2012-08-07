#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: ./catxyz.sh filename"
  exit
fi

for file in $1*.xyz; do
  N=`head -n 1 $file`
  cat $file >> $1_all_$N.xyz
done
