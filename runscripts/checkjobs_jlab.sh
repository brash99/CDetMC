#!/bin/sh
echo 'Running jobs:'
squeue -u brash | grep R | grep farm | wc -l
echo 'Pending jobs:'
squeue -u brash | grep PD | wc -l
