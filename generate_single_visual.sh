#!/bin/sh

cd batch
root -b -q 'GenParticles.C(13,50,1,0.0,0.0)'
cp AnaBarMC_Gen_1.root ../data/Gen_test1.root
cd ..
AnaBarMC macros/vis.mac
