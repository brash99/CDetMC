#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "AddNoise.C"
#include "ChiSquareTest.C"

using namespace std;

void noise_removal_test() {

    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 1);
    ChiSquareTest(3007);

    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 2);
    ChiSquareTest(3007);
    
    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 3);
    ChiSquareTest(3007);
    
    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 4);
    ChiSquareTest(3007);
    
    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 5);
    ChiSquareTest(3007);
    
    gSystem->Exec("cp data/AnaBarMC_3007_orig.root data/AnaBarMC_3007.root");
    AddNoise(3007, 10);
    ChiSquareTest(3007);

    return;
}

