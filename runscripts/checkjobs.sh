#!/bin/sh
echo 'Running jobs:'
squeue -u brash | grep R | grep product | grep jlab | wc -l
echo 'Pending jobs:'
squeue -u brash | grep PD | wc -l
