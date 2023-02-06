#!/bin/bash

if [[ $# -ne 1 ]]; then
	echo "illegal number of parameters"
	exit 
fi

export runlist=$1
export homedir=$HOME
export OUTPUT_DIR=$homedir/CDetOptical/data

echo "Reading runs from $runlist"

#submit a job for each file in filelist
for run in `cat $runlist`
do
    echo "The next run is $run"
    export RUN_NUMBER=$run
#    qsub AnaBarNeutron.sh 
    ./AnaBarNeutronUbuntu.sh
    sleep 1
done

echo "All finished ... :)"
