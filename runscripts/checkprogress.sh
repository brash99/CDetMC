#!/bin/sh
echo 'Last 10 lines of most recent slurm output file:'
ls -altr ../batch/slurm-*.out | tail -1 | awk '{print $9}' | xargs tail
