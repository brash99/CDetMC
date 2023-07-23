#!/bin/bash

#Initialize values
fPDGVal=${1:-"13"}
start=${2:-7001}
end=${3:-7064}
nevents=${4:-100}
xvalue=${5:-0}
zvalue=${6:-0}

#The generation code
echo "Generating input ROOT files for given particle ..."
cd batch
for (( i=$start ; i<=$end ; i++ ))
do
root -b -q 'GenParticles.C('$fPDGVal','$nevents','$i','$xvalue','$zvalue')'
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
ssh jlabanalysis << EOF
cd ~/CDetMC/batch
./sendjobs.sh runlist_test
#echo "Done"
EOF





