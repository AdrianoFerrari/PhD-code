#!/bin/bash

for job in *.pbs; do
  qsub $job
done
