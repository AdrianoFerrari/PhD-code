#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: ./batcher.sh runmin batchmin"
  exit
fi

PER=`echo "$2/(1.15*$1)" | bc`

cat ompjob* | grep main > bodies
split -l $PER bodies split_

rm ompjob*
rm bodies

i=0
for file in split_*; do
  cat head $file > ompjob$i.pbs
  let i=i+1
done

rm split_*
