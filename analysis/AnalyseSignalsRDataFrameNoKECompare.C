#include <iostream>

using namespace std;
using RNode = ROOT::RDF::RNode;

std::vector<RNode> v;
int global_run_number1;
int global_run_number2;
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;

int MaxPMTNo = 50000;
int MaxPMTHits = 1000;
float Finger_Edep_Max = 10.0;
float AnaBar_Edep_Max = 10.0;
float pedastel_sigma = 2.9;
int Detector_Offset = 2560;
int Detector_PMT_Offset = 2500;
int AnaBar_Offset = 30000;
int AnaBar_PMT_Offset = 0;

int Finger_NPhotons_Max = 250;
int AnaBar_NPhotons_Max = 200;

const int NUMPADDLE = 14;
const int NUMBARS = 14;
const int NUMMODULES = 3;
const int NUMSIDES = 2;
const int NUMLAYERS = 2;

const int NDET = NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;

int NMaxPMT = 14;

bool getTrigger(int Detector_Nhits, int* Detector_id) {

    bool tophit = false;
    bool bottomhit = false;
    bool fhit = false;
    bool ahit = false;
    bool trigger = false;
    for (int j=0; j<Detector_Nhits; j++) {
        //std::cout << "Detector id = " << Detector_id[j] << std::endl;
        if ((Detector_id[j] == Detector_Offset || Detector_id[j] == Detector_Offset+1) && !tophit) {
            tophit = true;
            //std::cout << "Top hit" << Detector_id[j] << std::endl;
        }
        if ((Detector_id[j] == Detector_Offset+2 || Detector_id[j] == Detector_Offset+3) && !bottomhit) {
            bottomhit = true;
            //std::cout << "Bottom hit" << Detector_id[j] << std::endl;
        }
        for (int ibar=0; ibar<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; ibar++){
            if (Detector_id[j] == AnaBar_Offset + ibar) {
                ahit = true;
            }
        }
        //if (tophit && bottomhit) {
        if (ahit) {
            fhit = true;
            trigger = true;
        }
    }

    return trigger;
}

bool getTrigger2(bool trigger, float fNewTheta) {

    bool trigger2 = false;
    if (fNewTheta > Theta_min_cut) {
        trigger2 = true;
    }
    return trigger2;
}

bool getTrigger3(bool trigger, float fNewTheta) {

    bool trigger3 = false;
    if (fNewTheta > ThetaVerticalCut) {
        trigger3 = true;
    }
    return trigger3;
}

float getMass(int Prim_pdg) {

    float fMass;
    if (Prim_pdg == 11) {
        fMass = 0.511;
    } else {
        if (Prim_pdg == 13) {
            fMass = 105.7;
        } else {
            if (Prim_pdg == 2212) {
                fMass = 938.28;
            } else {
                fMass = 939.65;
            }
        }
    }
    return fMass;
}

float getMomentum(float Prim_E, float fMass) {
    return sqrt(Prim_E*Prim_E - fMass*fMass);
}

float getPx(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Sin(Prim_Th)*TMath::Cos(Prim_Ph);
}

float getPy(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Sin(Prim_Th)*TMath::Sin(Prim_Ph);
}

float getPz(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Cos(Prim_Th);
}

float getNewTheta(float fMomentum, float fPy) {
    return TMath::ACos(fPy/fMomentum);
}

float getNewPhi(float fMomentum, float fPx, float fPz) {
    float fNewPhi;
    if (fPx < 0.0) {
        fNewPhi = TMath::ATan(fPz/fPx) + TMath::Pi();
    } else {
        if (fPx > 0.0 && fPz < 0.0) {
            fNewPhi = TMath::ATan(fPz/fPx) + TMath::TwoPi();
        } else {
            fNewPhi = TMath::ATan(fPz/fPx);
        }
    }
    return fNewPhi;
}


int getAnaBarMult(bool trigger, int* PMT_Nphotons) {

    int imult = 0;
    float temp;
    TRandom3* fRand = new TRandom3(-1);

    for (int icount = 0;icount < NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
        temp = PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
        if (temp>3.0*pedastel_sigma) {
            imult++;
        }
    }

    return imult;
}

std::vector<float> getFingerXVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_x, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3) && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_x[j]);
            }
        }
    }
    return v;
}
std::vector<float> getFingerYVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_y, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3) && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_y[j]);
            }
        }
    }
    return v;
}
std::vector<float> getFingerZVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_z, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3) && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_z[j]);
            }
        }
    }
    return v;
}

std::vector<float> getFingerTVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_t, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3) && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_t[j]);
            }
        }
    }
    return v;
}

std::vector<float> getAnaBarXVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_x, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_x[j]);
            }
        }
    }
    return v;
}
std::vector<float> getAnaBarYVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_y, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_y[j]);
            }
        }
    }
    return v;
}

std::vector<float> getAnaBarZVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_z, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_z[j]);
            }
        }
    }
    return v;
}
std::vector<float> getAnaBarTVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_t, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_t[j]);
            }
        }
    }
    return v;
}

std::vector<int> getFingerID(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {

    std::vector<int> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3)) {
                v.push_back(Detector_id[j]);
            }
        }
    }
    return v;
}

std::vector<int> getFingerPDG(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {

    std::vector<int> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if ((Detector_id[j] >= Detector_Offset && Detector_id[j] <= Detector_Offset+3)) {
                v.push_back(Detector_pdg[j]);
            }
        }
    }
    return v;
}

std::vector<int> getAnaBarID(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {

    std::vector<int> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS) {
                v.push_back(Detector_id[j]);
            }
        }
    }
    return v;
}

std::vector<int> getAnaBarPDG(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {

    std::vector<int> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS) {
                v.push_back(Detector_pdg[j]);
            }
        }
    }
    return v;
}

std::vector<float> getFingerPMTNPhotons(bool trigger, int* PMT_Nphotons) {

    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmt0tot = 0;


    if (trigger) {
        for (Int_t icount = Detector_PMT_Offset;icount<Detector_PMT_Offset+2;icount++){
            //std::cout << "getFingerPMTNphotons: " << icount << " " << PMT_Nphotons[icount] << std::endl;
            pmt0tot += PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
        }
    }

    v.push_back(pmt0tot);

    return v;
}

std::vector<float> getAnaBarPMTNPhotons(bool trigger, int* PMT_Nphotons) {

    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttot[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];

    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>0) {
                //std::cout << "getFingerPMTNphotons: " << icount << " " << PMT_Nphotons[icount] << std::endl;
                pmttot[icount] = PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
                v.push_back(pmttot[icount]);
            }
        }
    }
    return v;
}

std::vector<float> getFingerPMTID(bool trigger, int* PMT_Nphotons) {

    std::vector<float> v;

    if (trigger) {
        for (Int_t icount = Detector_PMT_Offset;icount<Detector_PMT_Offset+2;icount++){
            if (PMT_Nphotons[icount]>0) {
                //std::cout << "getFingerPMTID: " << icount << " " << PMT_Nphotons[icount] << std::endl;
                v.push_back(icount);
            }
        }
    }

    return v;
}

std::vector<float> getAnaBarPMTID(bool trigger, int* PMT_Nphotons) {

    std::vector<float> v;

    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>0) {
                //std::cout << "getAnaBarPMTID: " << icount << " " << PMT_Nphotons[icount] << std::endl;
                v.push_back(icount);
            }
        }
    }

    return v;
}

std::vector<float> getAnaBarNPhotonsTotal(bool trigger, int* PMT_Nphotons) {

    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttot = 0;

    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>0) {
                pmttot = pmttot + PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
            }
        }
    }

    v.push_back(pmttot);

    return v;
}

std::vector<float> getFingerEd(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {

    std::vector<float> v;
    float edep0tot = 0;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset || Detector_id[j] == Detector_Offset+1) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edep0tot += Detector_Ed[j];
                }
            }
        }
    }

    v.push_back(edep0tot);

    return v;
}

std::vector<float> getAnaBarEd(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {

    std::vector<float> v;
    float edeptot[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];

    for (int j=0; j<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; j++) {
        edeptot[j]=0.0;
    }

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_Ed[j] > 0.0 && Detector_id[j] >= AnaBar_Offset && Detector_id[j] <= AnaBar_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edeptot[Detector_id[j]-AnaBar_Offset] += Detector_Ed[j];
                    //std::cout << "j: " << j << "  index: " << Detector_id[j]-AnaBar_Offset << " energy: " << Detector_Ed[j] << std::endl;
                } else {
                    if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
                        edeptot[Detector_id[j]-AnaBar_Offset] += Detector_Ed[j];
                    }
                }
            }
        }
    }

    for (int j=0; j<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; j++) {
        if (edeptot[j] > 0.0) {
            v.push_back(edeptot[j]);
        }
    }

    return v;
}

std::vector<float> getAnaBarEdTotal(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {

    std::vector<float> v;
    float edeptotal = 0;
    float edeptot[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];

    for (int j=0; j<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; j++) {
        edeptot[j]=0.0;
    }

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_Ed[j] > 0.0 && Detector_id[j] >= AnaBar_Offset  && Detector_id[j] <= AnaBar_Offset + NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edeptot[Detector_id[j]-AnaBar_Offset] += Detector_Ed[j];
                } else {
                    if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
                        edeptot[Detector_id[j]-AnaBar_Offset] += Detector_Ed[j];
                    }
                }
            }
        }
    }

    for (int j=0; j<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; j++) {
        if (edeptot[j] > 0.0) {
            edeptotal += edeptot[j];
        }
    }

    v.push_back(edeptotal);

    return v;
}

//std::vector<RNode> AnalyseSignalsRDataFrameNoKECompare(int run_number1 = 4000, int run_number2 = 4001) {
void AnalyseSignalsRDataFrameNoKECompare(int run_number1 = 4000, int run_number2 = 4001) {

        //std::vector<RNode> v;
	global_run_number1 = run_number1;
	global_run_number2 = run_number2;
	std::cout << run_number1 << std::endl;
	std::cout << run_number2 << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);
	auto fileName1 = "data/AnaBarMC_"+std::to_string(run_number1)+".root";
	auto fileName2 = "data/AnaBarMC_"+std::to_string(run_number2)+".root";
        std::cout << fileName1 << std::endl;
        std::cout << fileName2 << std::endl;
	auto treeName = "T";

	ROOT::RDataFrame d1(treeName,fileName1);
	ROOT::RDataFrame d2(treeName,fileName2);

	//auto entries = d.Count();
	//cout << *entries << " entries in Tree with no filter" << endl;

	auto fdf1 = d1.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])")
       			.Define("fMass", "getMass(Prim_pdg)")
       			.Define("fMomentum","getMomentum(Prim_E,fMass)")
       			.Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fNewTheta", "getNewTheta(fMomentum,fPy)")
       			.Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)")
       			.Define("trigger2", "getTrigger2(trigger,fNewTheta)")
       			.Define("trigger3", "getTrigger3(trigger,fNewTheta)")
       			.Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       	                .Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
                        .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");

	auto fdf2 = d2.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])")
       			.Define("fMass", "getMass(Prim_pdg)")
       			.Define("fMomentum","getMomentum(Prim_E,fMass)")
       			.Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fNewTheta", "getNewTheta(fMomentum,fPy)")
       			.Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)")
       			.Define("trigger2", "getTrigger2(trigger,fNewTheta)")
       			.Define("trigger3", "getTrigger3(trigger,fNewTheta)")
       			.Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
                        .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");
	
        //auto entries2 = fdf.Count();
	//cout << *entries2 << " entries in Expanded Dataframe with no filter" << endl;

	auto triggers1 = fdf1.Filter("trigger==true").Count();
	cout << *triggers1<< " File 1 entries passed Main trigger" << endl;
	auto triggers2 = fdf2.Filter("trigger==true").Count();
	cout << *triggers2<< " File 2 entries passed Main trigger" << endl;

	auto fdft1 = fdf1.Filter("trigger==true");
	auto fdft2 = fdf2.Filter("trigger==true");
        v.push_back(fdft1);
        v.push_back(fdft2);


	//return v;
}

TCanvas* plotC1(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

  auto hFingerX1 = v[0].Histo1D("fingerXVec");
  auto hFingerY1 = v[0].Histo1D("fingerYVec");
  auto hFingerZ1 = v[0].Histo1D("fingerZVec");
  auto hFingerT1 = v[0].Histo1D("fingerTVec");
  auto hFingerX2 = v[1].Histo1D("fingerXVec");
  auto hFingerY2 = v[1].Histo1D("fingerYVec");
  auto hFingerZ2 = v[1].Histo1D("fingerZVec");
  auto hFingerT2 = v[1].Histo1D("fingerTVec");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);

  c1->cd(1);
  hFingerX1->Draw();
  c1->cd(2);
  hFingerY1->Draw();
  c1->cd(3);
  hFingerZ1->Draw();
  c1->cd(4);
  hFingerT1->Draw();
  c1->cd(1);
  hFingerX2->Draw("SAME");
  c1->cd(2);
  hFingerY2->Draw("SAME");
  c1->cd(3);
  hFingerZ2->Draw("SAME");
  c1->cd(4);
  hFingerT2->Draw("SAME");

  c1->DrawClone();
  c1->Print("plots/c1.pdf");

  return c1;

}

TCanvas* plotC2(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

        auto hPrimE1 = v[0].Histo1D({"h1", "E", 100, 0.0, 7000.0},"Prim_E");
        auto hPrimTh1 = v[0].Histo1D({"h1", "TH", 100, -0.1, 3.2},"fNewTheta");
        auto hPrimPh1 = v[0].Histo1D({"h1", "PH", 100, -0.1, 6.4},"fNewPhi");
        auto hPrimPdg1 = v[0].Histo1D({"h1", "PDG", 100, 0.0, 20.0},"Prim_pdg");
        auto hPrimE2 = v[1].Histo1D({"h1", "E", 100, 0.0, 7000.0},"Prim_E");
        auto hPrimTh2 = v[1].Histo1D({"h1", "TH", 100, -0.1, 3.2},"fNewTheta");
        auto hPrimPh2 = v[1].Histo1D({"h1", "PH", 100, -0.1, 6.4},"fNewPhi");
        auto hPrimPdg2 = v[1].Histo1D({"h1", "PDG", 100, 0.0, 20.0},"Prim_pdg");

	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(2,2,0.01,0.01,0);

	c2->cd(1);
	hPrimE2->Draw();
	c2->cd(2);
	hPrimTh2->Draw();
	c2->cd(3);
	hPrimPh2->Draw();
	c2->cd(4);
	hPrimPdg2->Draw();
	c2->cd(1);
	hPrimE1->Draw("SAME");
	c2->cd(2);
	hPrimTh1->Draw("SAME");
	c2->cd(3);
	hPrimPh1->Draw("SAME");
	c2->cd(4);
	hPrimPdg1->Draw("SAME");

	c2->DrawClone();
	c2->Print("plots/c2RA.pdf");

	return c2;

}


TCanvas* plotC3(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hDetectorNhits1 = v[0].Histo1D("Detector_Nhits");
	auto hDetectorPdg1 = v[0].Histo1D("anaBarPDG");
	auto hDetectorID1 = v[0].Histo1D("anaBarID");
	auto hFingerPdg1 = v[0].Histo1D("fingerPDG");
	auto hFingerID1 = v[0].Histo1D("fingerID");
	auto hPMTID1 = v[0].Histo1D("PMT_id");
	auto hAnaBarPMTID1 = v[0].Histo1D("anaBarPMTID");
	auto hFingerPMTID1 = v[0].Histo1D("fingerPMTID");
	
        auto hDetectorNhits2 = v[1].Histo1D("Detector_Nhits");
	auto hDetectorPdg2 = v[1].Histo1D("anaBarPDG");
	auto hDetectorID2 = v[1].Histo1D("anaBarID");
	auto hFingerPdg2 = v[1].Histo1D("fingerPDG");
	auto hFingerID2 = v[1].Histo1D("fingerID");
	auto hPMTID2 = v[1].Histo1D("PMT_id");
	auto hAnaBarPMTID2 = v[1].Histo1D("anaBarPMTID");
	auto hFingerPMTID2 = v[1].Histo1D("fingerPMTID");

	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	c3->Divide(3,3,0.01,0.01,0);

	c3->cd(1);
	hDetectorNhits1->Draw();
	c3->cd(2);
	hFingerPdg1->Draw();
	c3->cd(3);
	hDetectorPdg1->Draw();
	c3->cd(4);
	hFingerID1->Draw();
	c3->cd(5);
	hDetectorID1->Draw();
	c3->cd(6);
	hPMTID1->Draw();
	c3->cd(7);
	hFingerPMTID1->Draw();
	c3->cd(8);
	hAnaBarPMTID1->Draw();
	c3->cd(1);
	hDetectorNhits2->Draw("SAME");
	c3->cd(2);
	hFingerPdg2->Draw("SAME");
	c3->cd(3);
	hDetectorPdg2->Draw("SAME");
	c3->cd(4);
	hFingerID2->Draw("SAME");
	c3->cd(5);
	hDetectorID2->Draw("SAME");
	c3->cd(6);
	hPMTID2->Draw("SAME");
	c3->cd(7);
	hFingerPMTID2->Draw("SAME");
	c3->cd(8);
	hAnaBarPMTID2->Draw("SAME");

	c3->DrawClone();
	c3->Print("plots/c3RA.pdf");

	return c3;

}


TCanvas* plotC4(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hFingerEd1 = v[0].Histo1D({"h1", "Finger EDep", 100, 0.0, 10.0},"fingerEd");
	auto hFingerPMTNphot1 = v[0].Histo1D({"h1", "Finger Npe", 100, 0.0, 20.0},"fingerPMTNPhotons");
	auto hAnaBarPMTNphot1 = v[0].Histo1D({"h1", "AnaBar Npe", 100, 0.0, 200.0},"anaBarPMTNPhotons");
	auto hAnaBarEd1 = v[0].Histo1D({"h1", "AnaBar EDep", 100, 0.0, 10.0},"anaBarEd");
	auto hFingerEd2 = v[1].Histo1D({"h1", "Finger EDep", 100, 0.0, 10.0},"fingerEd");
	auto hFingerPMTNphot2 = v[1].Histo1D({"h1", "Finger Npe", 100, 0.0, 20.0},"fingerPMTNPhotons");
	auto hAnaBarPMTNphot2 = v[1].Histo1D({"h1", "AnaBar Npe", 100, 0.0, 200.0},"anaBarPMTNPhotons");
	auto hAnaBarEd2 = v[1].Histo1D({"h1", "AnaBar EDep", 100, 0.0, 10.0},"anaBarEd");

	TCanvas *c4 = new TCanvas("c4","c4",800,800);

	c4->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.50,0.99);
	pad1->Draw();
	pad1->cd();
	hFingerEd2->GetXaxis()->SetRangeUser(1.0,10);
	hFingerEd2->Draw();
	hFingerEd1->GetXaxis()->SetRangeUser(1.0,10);
	hFingerEd1->Draw("SAME");

	c4->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
	pad2->Draw();
	pad2->cd();
	hFingerPMTNphot2->GetXaxis()->SetRangeUser(-10,250);
	hFingerPMTNphot2->Draw();
	hFingerPMTNphot1->GetXaxis()->SetRangeUser(-10,250);
	hFingerPMTNphot1->Draw("SAME");

	c4->cd();
	TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.50,0.50);
	//pad3->SetLogy();
	pad3->Draw();
	pad3->cd();
	hAnaBarEd2->GetXaxis()->SetRangeUser(1.0,10);
	hAnaBarEd2->Draw();
	hAnaBarEd1->GetXaxis()->SetRangeUser(1.0,10);
	hAnaBarEd1->Draw("SAME");

	c4->cd();
	TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.50);
	//pad4->SetLogy();
	pad4->Draw();
	pad4->cd();
	hAnaBarPMTNphot2->GetXaxis()->SetRangeUser(-20,180);
	hAnaBarPMTNphot2->Draw();
	hAnaBarPMTNphot1->GetXaxis()->SetRangeUser(-20,180);
	hAnaBarPMTNphot1->Draw("SAME");

	c4->DrawClone();
	c4->Print("plots/c4RA.pdf");

	return c4;
}


TCanvas* plotC5(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hAnaBarX1 = v[0].Histo1D({"h1", "AnaBar X", 100, -800.0, 800.0},"anaBarXVec");
	auto hAnaBarY1 = v[0].Histo1D({"h1", "AnaBar Y", 100, -400.0, 100.0},"anaBarYVec");
	auto hAnaBarZ1 = v[0].Histo1D({"h1", "AnaBar Z", 100, -2000.0, 2000.0},"anaBarZVec");
	auto hAnaBarT1 = v[0].Histo1D({"h1", "AnaBar T", 100, 0.0, 10.0},"anaBarTVec");
	auto hAnaBarX2 = v[1].Histo1D({"h1", "AnaBar X", 100, -800.0, 800.0},"anaBarXVec");
	auto hAnaBarY2 = v[1].Histo1D({"h1", "AnaBar Y", 100, -400.0, 100.0},"anaBarYVec");
	auto hAnaBarZ2 = v[1].Histo1D({"h1", "AnaBar Z", 100, -2000.0, 2000.0},"anaBarZVec");
	auto hAnaBarT2 = v[1].Histo1D({"h1", "AnaBar T", 100, 0.0, 10.0},"anaBarTVec");

	TCanvas *c5 = new TCanvas("c5","c5",800,800);
	c5->Divide(2,2,0.01,0.01,0);

	c5->cd(1);
	hAnaBarX1->Draw("SAME");
	c5->cd(2);
	hAnaBarY1->Draw("SAME");
	c5->cd(3);
	hAnaBarZ1->Draw("SAME");
	c5->cd(4);
	hAnaBarT1->Draw("SAME");
	c5->cd(1);
	hAnaBarX2->Draw("SAME");
	c5->cd(2);
	hAnaBarY2->Draw("SAME");
	c5->cd(3);
	hAnaBarZ2->Draw("SAME");
	c5->cd(4);
	hAnaBarT2->Draw("SAME");

	c5->DrawClone();
	c5->Print("plots/c5RA.pdf");

	return c5;

}

TCanvas* plotC6(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hE1vsE21 = v[0].Histo2D({"h2", "E1 vs E2", 100, 0.01, 3.0, 100, 0.01, 30.0},"fingerEd","anaBarEdTotal");
	auto hE1vsE22 = v[1].Histo2D({"h2", "E1 vs E2", 100, 0.01, 3.0, 100, 0.01, 30.0},"fingerEd","anaBarEdTotal");

	TCanvas *c6 = new TCanvas("c6","c6",800,800);
	c6->Divide(2,1,0.01,0.01,0);

	c6->cd(1);
	hE1vsE21->Draw("COLZ");
	c6->cd(2);
	hE1vsE22->Draw("COLZ");

	c6->DrawClone();
	c6->Print("plots/c6RA.pdf");

	return c6;

}

TCanvas* plotC7(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hFinger_Edep_vs_Nphot1 = v[0].Filter("trigger2").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot1 = v[0].Filter("trigger2").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1 = v[0].Filter("trigger2").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");
	auto hFinger_Edep_vs_Nphot2 = v[1].Filter("trigger2").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot2 = v[1].Filter("trigger2").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot2 = v[1].Filter("trigger2").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c7 = new TCanvas("c7","c7",800,800);
	c7->Divide(4,2,0.01,0.01,0);

	c7->cd(1);
	hFinger_Edep_vs_Nphot1->Draw("COLZ");
	c7->cd(2);
	hAnaBar_Edep_vs_Nphot1->Draw("COLZ");
	c7->cd(3);
	hNphot0_vs_Nphot1->Draw("COLZ");
	c7->cd(4);
	TProfile *prof1 = hAnaBar_Edep_vs_Nphot1->ProfileX();
	prof1->Fit("pol1");
	c7->cd(5);
	hFinger_Edep_vs_Nphot2->Draw("COLZ");
	c7->cd(6);
	hAnaBar_Edep_vs_Nphot2->Draw("COLZ");
	c7->cd(7);
	hNphot0_vs_Nphot2->Draw("COLZ");
	c7->cd(8);
	TProfile *prof2 = hAnaBar_Edep_vs_Nphot2->ProfileX();
	prof2->Fit("pol1");

	c7->DrawClone();
	c7->Print("plots/c7RA.pdf");

	return c7;

}

TCanvas* plotC8(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hFinger_Edep_vs_Nphot1 = v[0].Filter("trigger3").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot1 = v[0].Filter("trigger3").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1 = v[0].Filter("trigger3").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");
	auto hFinger_Edep_vs_Nphot2 = v[1].Filter("trigger3").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot2 = v[1].Filter("trigger3").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot2 = v[1].Filter("trigger3").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c8 = new TCanvas("c8","c8",800,800);
	c8->Divide(4,2,0.01,0.01,0);

	c8->cd(1);
	hFinger_Edep_vs_Nphot1->Draw("COLZ");
	c8->cd(2);
	hAnaBar_Edep_vs_Nphot1->Draw("COLZ");
	c8->cd(3);
	hNphot0_vs_Nphot1->Draw("COLZ");
	c8->cd(4);
	TProfile *prof1 = hAnaBar_Edep_vs_Nphot1->ProfileX();
	prof1->Fit("pol1");
	c8->cd(5);
	hFinger_Edep_vs_Nphot2->Draw("COLZ");
	c8->cd(6);
	hAnaBar_Edep_vs_Nphot2->Draw("COLZ");
	c8->cd(7);
	hNphot0_vs_Nphot2->Draw("COLZ");
	c8->cd(8);
	TProfile *prof2 = hAnaBar_Edep_vs_Nphot2->ProfileX();
	prof2->Fit("pol1");

	c8->DrawClone();
	c8->Print("plots/c8RA.pdf");

	return c8;

}

TCanvas* plotC11(){
    
    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

    auto hAnaBarMult1 = v[0].Histo1D({"h1", "CDet Multiplicity", 20, 0, 20},"imult");
    auto hAnaBarMult2 = v[1].Histo1D({"h1", "CDet Multiplicity", 20, 0, 20},"imult");

    TCanvas* c11 = new TCanvas("c11", "c11", 800,800);
    c11->Divide(1,1, 0.01, 0.01, 0);

    c11->cd(1);
    hAnaBarMult1->Draw();
    hAnaBarMult2->Draw("SAME");

    c11->DrawClone();
    c11->Print("plots/c11RA.pdf");

    return c11;

}

TCanvas* plotC12(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hPrimPx1 = v[0].Histo1D({"h1", "Px", 100, -1000, 1000},"fPx");
	auto hPrimPy1 = v[0].Histo1D({"h1", "Px", 100, -8000, 8000},"fPy");
	auto hPrimPz1 = v[0].Histo1D({"h1", "Px", 100, -2000, 2000},"fPz");
	auto hPrimPx2 = v[1].Histo1D({"h1", "Px", 100, -1000, 1000},"fPx");
	auto hPrimPy2 = v[1].Histo1D({"h1", "Px", 100, -8000, 8000},"fPy");
	auto hPrimPz2 = v[1].Histo1D({"h1", "Px", 100, -2000, 2000},"fPz");

	TCanvas *c12 = new TCanvas("c12","c12",800,800);
	c12->Divide(2,2,0.01,0.01,0);

	c12->cd(1);
	hPrimPx1->Draw("SAME");
	c12->cd(2);
	hPrimPy1->Draw("SAME");
	c12->cd(3);
	hPrimPz1->Draw("SAME");
	c12->cd(1);
	hPrimPx2->Draw("SAME");
	c12->cd(2);
	hPrimPy2->Draw("SAME");
	c12->cd(3);
	hPrimPz2->Draw("SAME");

	c12->DrawClone();
	c12->Print("plots/c12.pdf");

	return c12;

}


TCanvas* plotC13(){

    //std::vector<RNode> v = AnalyseSignalsRDataFrameNoKECompare(global_run_number1,global_run_number2);

	auto hPx_vs_x1 = v[0].Filter("trigger2").Histo2D({"h33", "Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0},"anaBarXVec","fPx");
	auto hPz_vs_z1 = v[0].Filter("trigger2").Histo2D({"h34", "Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0},"anaBarZVec","fPz");
	auto hz_vs_x1 = v[0].Filter("trigger2").Histo2D({"h35", "z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0},"anaBarXVec","anaBarZVec");
	auto hPx_vs_x2 = v[1].Filter("trigger2").Histo2D({"h33", "Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0},"anaBarXVec","fPx");
	auto hPz_vs_z2 = v[1].Filter("trigger2").Histo2D({"h34", "Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0},"anaBarZVec","fPz");
	auto hz_vs_x2 = v[1].Filter("trigger2").Histo2D({"h35", "z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0},"anaBarXVec","anaBarZVec");

	TCanvas *c13 = new TCanvas("c13","c13",800,800);
	c13->Divide(3,2,0.01,0.01,0);

	c13->cd(1);
	hPx_vs_x1->Draw("COLZ");
	c13->cd(2);
	hPz_vs_z1->Draw("COLZ");
	c13->cd(3);
	hz_vs_x1->Draw("COLZ");
	c13->cd(4);
	hPx_vs_x2->Draw("COLZ");
	c13->cd(5);
	hPz_vs_z2->Draw("COLZ");
	c13->cd(6);
	hz_vs_x2->Draw("COLZ");

	c13->DrawClone();
	c13->Print("plots/c13.pdf");

	return c13;

}

TCanvas* plotC14(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimX1 = v[0].Histo1D({"h99","X_vtx", 100, -80,80},"Prim_X");
	auto hPrimY1 = v[0].Histo1D({"h99","Y_vtx", 100, 0,40},"Prim_Y");
	auto hPrimZ1 = v[0].Histo1D({"h99","Z_vtx", 100, -200,200},"Prim_Z");
	auto hPrimXZ1 = v[0].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");
	auto hPrimX2 = v[1].Histo1D({"h99","X_vtx", 100, -80,80},"Prim_X");
	auto hPrimY2 = v[1].Histo1D({"h99","Y_vtx", 100, 0,40},"Prim_Y");
	auto hPrimZ2 = v[1].Histo1D({"h99","Z_vtx", 100, -200,200},"Prim_Z");
	auto hPrimXZ2 = v[1].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");
	
	TCanvas *c14 = new TCanvas("c14","c14",800,800);
	c14->Divide(3,2,0.01,0.01,0);

	c14->cd(1);
	hPrimX1->Draw();
	c14->cd(2);
	hPrimY1->Draw();
	c14->cd(3);
	hPrimZ1->Draw();
	c14->cd(4);
	hPrimXZ1->Draw("COLZ");
	c14->cd(1);
	hPrimX2->Draw("SAME");
	c14->cd(2);
	hPrimY2->Draw("SAME");
	c14->cd(3);
	hPrimZ2->Draw("SAME");
	c14->cd(5);
	hPrimXZ2->Draw("COLZ");

	c14->DrawClone();
	c14->Print("plots/c14.pdf");

	return c14;

}
