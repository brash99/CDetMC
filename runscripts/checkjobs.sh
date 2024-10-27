#!/bin/sh
echo 'Running jobs:'
squeue -u $USER | grep R | grep product | grep jlab | wc -l
echo 'Pending jobs:'
squeue -u $USER | grep PD | wc -l
