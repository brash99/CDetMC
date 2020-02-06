# AnaBarMC
To compile the code type:
make or gmake

To run the event generator type:
root -q GenCosmics.C++

This creates a root file in the data/ directory

To run the simulation interactively type:

AnaBarMC [return]
then at the prompt type:
/control/execute macros/vis.mac

To run the simaultion in batch mode type:

AnaBarMC macros/batch.mac

---Tanner's Comments---

To Run Simulation

>Go into ~/CDetOptical/batch directory and type: "root generate_80.C" to create the primary particle rootfiles that are then store in ~/CDetOptical/batch/data
                                                                                  
>Type: "./sendjobs.sh runlist80" to send the primary particle rootfiles that are   listed in the according runlist file that is sent as an argumant. The output data will be stored in ~/CDetOptical/data
                                                                                  
Changing the Primary Particle

>Switch the Generator Function in ~/CDetOptical/batch/GenCosmics.C at line 100

>Change PrimaryParticleID in ~/CDetOptical/AnalyseSignals.C at line 28

Detector Position and Coordinate Allignment

>X: Length of detector (~50cm for bar(y= to y=) and ~10cm for trigger paddle(y=-5cm to y=5cm)) 

>Y: Height/Thickness of detector (particles are generated @y=5cm, the Trigger Paddle exitsts from y=4cm to y=3cm, and the bar exists from y=2cm to y=-2cm) 

>Z: Width of detector (depends on #bars, 14 paddles per bar @ 0.5cm per bar) Goes from z=0 to z=-8cm for one bar (14 paddles @ 0.5cm per paddle plus some extra for the thicckness of mylar sheets between each paddle) 

Looking at Data via ROOT

>Open TBrowser and look at TTree sets

>Run AnalyseSignals('run#') and then call plotC1() -> plotC12() 
