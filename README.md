# CDetMC
To compile the code type:
make or gmake

To run the event generator type:

root -b -q 'batch/GenParticles.C(13,100,7000,0.0,0.0)'

where the arguments are:
ParticleType, Number of Events, Run Number, x position, z position

This creates a root file in the data/ directory

To run the simulation interactively type:

AnaBarMC [return]
then at the prompt type:
/control/execute macros/vis_single.mac

To run the simulation in batch mode type:

./Automate.sh 13 7000 7010 100 0.0 0.0

where the arguments are:
ParticleType, Run Number Start, Run Number End, Number of Events per run, x position, z position

Note:  
For running on stand-alone Ubuntu systems, there are corresponding scripts with *Ubuntu* extensions.
For running on stand-alone Mac systems, there are corresponding scripts with *Mac* extensions.

