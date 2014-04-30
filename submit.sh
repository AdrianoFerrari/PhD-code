#!/bin/bash

for job in *.pbs; do
  msub -q sw $job
done
