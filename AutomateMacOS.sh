#!/bin/bash

#Initialize values
fPDGVal=${1:-"13"}
start=${2:-7001}
end=${3:-7080}
nevents=${4:-96}

#The generation code
echo "Generating input ROOT files for given particle ..."
cd batch
for (( i=$start ; i<=$end ; i++ ))
do
root -b -q 'GenParticlesMacOS.C('$fPDGVal','$nevents','$i')'
done
cd ..
echo "Done"

#Generate the runlist
echo "Generating run list ..."
echo $start > batch/runlist_test
for (( i=$(($start + 1)) ; i<=$end ; i++ ))
do
echo "$i" >> batch/runlist_test
done
echo "Done"

#Simulation
echo "Submitting simulation jobs ..."
cd ~/CDetOptical/batch
./sendjobsmacos.sh runlist_test
echo "Done"





