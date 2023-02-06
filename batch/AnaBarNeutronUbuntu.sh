#!/bin/bash
export homedir=$HOME

export nevents=100
export tempdir=$homedir/CDetOptical/batch

export MACRO_PATH=$homedir/CDetOptical/macros/
export MCMACRO=$tempdir/AnaBarMC_$RUN_NUMBER.mac

echo "/control/macroPath $MACRO_PATH"	 	                         >   $MCMACRO
echo "/AnaBarMC/physics/addPhysics standard_opt3"                        >>   $MCMACRO
echo "/AnaBarMC/physics/optical 1"	                                 >>  $MCMACRO
echo "/AnaBarMC/detector/AnaBarXpos 0.00"	                         >>  $MCMACRO
echo "/run/initialize"                                                   >>  $MCMACRO
echo "/AnaBarMC/generator/Mode 1"                              >>  $MCMACRO
echo "/AnaBarMC/generator/InputFile $tempdir/data/AnaBarMC_Gen_$RUN_NUMBER.root" >>  $MCMACRO
echo "/AnaBarMC/analysis/setOutputFile $tempdir/rootfiles/AnaBarMC_$RUN_NUMBER.root" >>  $MCMACRO

cd $tempdir

export ROOTSYS=/cern/root/pro
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
AnaBarMC $MCMACRO #>& /dev/null
echo "****************** AnaBarMC Finished"

cp    ${tempdir}/rootfiles/"AnaBarMC_$RUN_NUMBER.root"   ${OUTPUT_DIR}/
rm -f ${tempdir}/rootfiles/"AnaBarMC_$RUN_NUMBER.root"
rm -f $MCMACRO
