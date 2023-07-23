#!/bin/bash

#Initialize values
fPDGVal=${1:-"13"}
start=${2:-7001}
nevents=${4:-96}

#Generate the runlist
echo "Generating run list ..."
echo $start > batch/runlist_test_single
echo "Done"

#Simulation
echo "Submitting simulation jobs ..."
cd ~/CDetMC/batch
./sendjobsangelo.sh runlist_test_single
echo "Done"





