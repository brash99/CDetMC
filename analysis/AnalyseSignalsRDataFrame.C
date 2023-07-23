#include <iostream>

using namespace std;
using RNode = ROOT::RDF::RNode;

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
        if (tophit && bottomhit) {
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

std::vector<float> getFingerPMTKE(bool trigger, int* PMT_Nphotons, float* PMT_KineticEnergy) {

    std::vector<float> v;

    for (Int_t iq=0; iq<1; iq++) {
    	for (Int_t jq=0; jq<1; jq++){
        	if(PMT_KineticEnergy[(Detector_PMT_Offset+iq)*MaxPMTHits+jq] > 0.001) {
            		v.push_back(1240.0/PMT_KineticEnergy[(Detector_PMT_Offset+iq)*MaxPMTHits+jq]);
        	};
	}
    }

    return v;
}

std::vector<float> getAnaBarPMTKE(bool trigger, int* PMT_Nphotons, float* PMT_KineticEnergy) {


    // A note to future generations:
    // the way that two-dimentional arrays get stored in the DataFrame is that they get converted to one-dimenisional arrays
    // 
    // PMT_KineticEnergy is a 2D array of size [MaxPMTNo][MaxPMTHits] = [50000][1000] in this version of the simulation.
    //
    // This gets converted to a 1D array of size [MaxPMTNo*MaxPMTHits].  So, for example, for PMT with detector id = 1234, for hit number 567, 
    // that correponds to PMT_KineticEnergy[1234*MaxPMTHits+567]
    
    std::vector<float> v;

    //for (Int_t iq=0; iq<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; iq++) {
    for (Int_t iq=0; iq<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES; iq++) { //for now, just look in the TOP LAYER to speed things up a bit
        for (Int_t jq=0; jq<1; jq++){ //just look at the FIRST hit on that PMT
            if(PMT_KineticEnergy[(AnaBar_PMT_Offset+iq)*MaxPMTHits+jq] > 0.001) {
                v.push_back(1240.0/PMT_KineticEnergy[(AnaBar_PMT_Offset+iq)*MaxPMTHits+jq]);
            };
        }
    }

    return v;
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

RNode AnalyseSignalsRDataFrame() {

	auto fileName = "data/AnaBarMC_7001.root";
	auto treeName = "T";

	ROOT::RDataFrame d(treeName,fileName);

	//auto entries = d.Count();
	//cout << *entries << " entries in Tree with no filter" << endl;

	auto fdf = d.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])")
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
       			.Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
			.Define("fingerPMTKE","getFingerPMTKE(trigger,&PMT_Nphotons[0],&PMT_KineticEnergy[0])")
       			.Define("anaBarPMTKE","getAnaBarPMTKE(trigger,&PMT_Nphotons[0],&PMT_KineticEnergy[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");

	//auto entries2 = fdf.Count();
	//cout << *entries2 << " entries in Expanded Dataframe with no filter" << endl;

	auto triggers = fdf.Filter("trigger==true").Count();
	cout << *triggers << " entries passed Main trigger" << endl;

	auto fdft = fdf.Filter("trigger==true");

	return fdft;
}

TCanvas* plotC1(){

  RNode fdft = AnalyseSignalsRDataFrame();

  auto hFingerX = fdft.Histo1D("fingerXVec");
  auto hFingerY = fdft.Histo1D("fingerYVec");
  auto hFingerZ = fdft.Histo1D("fingerZVec");
  auto hFingerT = fdft.Histo1D("fingerTVec");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);

  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();

  c1->DrawClone();
  c1->Print("plots/c1.pdf");

  return c1;

}

TCanvas* plotC2(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hPrimE = fdft.Histo1D("Prim_E");
	auto hPrimTh = fdft.Histo1D("fNewTheta");
	auto hPrimPh = fdft.Histo1D("fNewPhi");
	auto hPrimPdg = fdft.Histo1D("Prim_pdg");

	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(2,2,0.01,0.01,0);

	c2->cd(1);
	hPrimE->Draw();
	c2->cd(2);
	hPrimTh->Draw();
	c2->cd(3);
	hPrimPh->Draw();
	c2->cd(4);
	hPrimPdg->Draw();

	c2->DrawClone();
	c2->Print("plots/c2RA.pdf");

	return c2;

}


TCanvas* plotC3(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hDetectorNhits = fdft.Histo1D("Detector_Nhits");
	auto hDetectorPdg = fdft.Histo1D("anaBarPDG");
	auto hDetectorID = fdft.Histo1D("anaBarID");
	auto hFingerPdg = fdft.Histo1D("fingerPDG");
	auto hFingerID = fdft.Histo1D("fingerID");
	auto hPMTID = fdft.Histo1D("PMT_id");

	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	c3->Divide(3,2,0.01,0.01,0);

	c3->cd(1);
	hDetectorNhits->Draw();
	c3->cd(2);
	hFingerPdg->Draw();
	c3->cd(3);
	hDetectorPdg->Draw();
	c3->cd(4);
	hFingerID->Draw();
	c3->cd(5);
	hDetectorID->Draw();
	c3->cd(6);
	hPMTID->Draw();

	c3->DrawClone();
	c3->Print("plots/c3RA.pdf");

	return c3;

}


TCanvas* plotC4(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hFingerEd = fdft.Histo1D("fingerEd");
	auto hFingerPMTNphot = fdft.Histo1D("fingerPMTNPhotons");
	auto hAnaBarPMTNphot = fdft.Histo1D("anaBarPMTNPhotons");
	auto hAnaBarEd = fdft.Histo1D("anaBarEd");

	TCanvas *c4 = new TCanvas("c4","c4",800,800);

	c4->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.50,0.99);
	pad1->Draw();
	pad1->cd();
	hFingerEd->GetXaxis()->SetRangeUser(1.0,10);
	hFingerEd->Draw();

	c4->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
	pad2->Draw();
	pad2->cd();
	hFingerPMTNphot->GetXaxis()->SetRangeUser(-10,250);
	hFingerPMTNphot->Draw();

	c4->cd();
	TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.50,0.50);
	//pad3->SetLogy();
	pad3->Draw();
	pad3->cd();
	hAnaBarEd->GetXaxis()->SetRangeUser(1.0,10);
	hAnaBarEd->Draw();

	c4->cd();
	TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.50);
	//pad4->SetLogy();
	pad4->Draw();
	pad4->cd();
	hAnaBarPMTNphot->GetXaxis()->SetRangeUser(-20,180);
	hAnaBarPMTNphot->Draw();

	c4->DrawClone();
	c4->Print("plots/c4RA.pdf");

	return c4;
}


TCanvas* plotC5(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hAnaBarX = fdft.Histo1D("anaBarXVec");
	auto hAnaBarY = fdft.Histo1D("anaBarYVec");
	auto hAnaBarZ = fdft.Histo1D("anaBarZVec");
	auto hAnaBarT = fdft.Histo1D("anaBarTVec");

	TCanvas *c5 = new TCanvas("c5","c5",800,800);
	c5->Divide(2,2,0.01,0.01,0);

	c5->cd(1);
	hAnaBarX->Draw();
	c5->cd(2);
	hAnaBarY->Draw();
	c5->cd(3);
	hAnaBarZ->Draw();
	c5->cd(4);
	hAnaBarT->Draw();

	c5->DrawClone();
	c5->Print("plots/c5RA.pdf");

	return c5;

}

TCanvas* plotC6(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hE1vsE2 = fdft.Histo2D({"h2", "E1 vs E2", 100, 0.01, 10.0, 100, 0.01, 30.0},"fingerEd","anaBarEdTotal");

	TCanvas *c6 = new TCanvas("c6","c6",800,800);
	c6->Divide(1,1,0.01,0.01,0);

	c6->cd(1);
	hE1vsE2->Draw("COLZ");

	c6->DrawClone();
	c6->Print("plots/c6RA.pdf");

	return c6;

}

TCanvas* plotC7(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hFinger_Edep_vs_Nphot = fdft.Filter("trigger2").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 250.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot = fdft.Filter("trigger2").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1 = fdft.Filter("trigger2").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 100.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c7 = new TCanvas("c7","c7",800,800);
	c7->Divide(2,2,0.01,0.01,0);

	c7->cd(1);
	hFinger_Edep_vs_Nphot->Draw("COLZ");
	c7->cd(2);
	hAnaBar_Edep_vs_Nphot->Draw("COLZ");
	c7->cd(3);
	hNphot0_vs_Nphot1->Draw("COLZ");
	c7->cd(4);
	TProfile *prof = hAnaBar_Edep_vs_Nphot->ProfileX();
	prof->Fit("pol1");

	c7->DrawClone();
	c7->Print("plots/c7RA.pdf");

	return c7;

}

TCanvas* plotC8(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hFinger_Edep_vs_NphotCut = fdft.Filter("trigger3").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 250.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_NphotCut = fdft.Filter("trigger3").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1Cut = fdft.Filter("trigger3").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 100.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c8 = new TCanvas("c8","c8",800,800);
	c8->Divide(2,2,0.01,0.01,0);

	c8->cd(1);
	hFinger_Edep_vs_NphotCut->Draw("COLZ");
	c8->cd(2);
	hAnaBar_Edep_vs_NphotCut->Draw("COLZ");
	c8->cd(3);
	hNphot0_vs_Nphot1Cut->Draw("COLZ");
	c8->cd(4);
	TProfile *profCut = hAnaBar_Edep_vs_NphotCut->ProfileX();
	profCut->Fit("pol1");

	c8->DrawClone();
	c8->Print("plots/c8RA.pdf");

	return c8;

}

TCanvas* plotC11(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hAnaBarMult = fdft.Histo1D("imult");

	TCanvas *c11 = new TCanvas("c11", "c11", 800,800);
	c11->Divide(1,1, 0.01, 0.01, 0);

	c11->cd(1);
	hAnaBarMult->Draw();

	c11->DrawClone();
	c11->Print("plots/c11RA.pdf");

	return c11;

}


TCanvas* plotC12(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hPrimPx = fdft.Histo1D("fPx");
	auto hPrimPy = fdft.Histo1D("fPy");
	auto hPrimPz = fdft.Histo1D("fPz");

	TCanvas *c12 = new TCanvas("c12","c12",800,800);
	c12->Divide(2,2,0.01,0.01,0);

	c12->cd(1);
	hPrimPx->Draw();
	c12->cd(2);
	hPrimPy->Draw();
	c12->cd(3);
	hPrimPz->Draw();

	c12->DrawClone();
	c12->Print("plots/c12.pdf");

	return c12;

}

TCanvas* plotC10(){

	RNode fdft = AnalyseSignalsRDataFrame();

	auto hFingerPMTKE = fdft.Histo1D("fingerPMTKE");
	auto hAnaBarPMTKE = fdft.Histo1D("anaBarPMTKE");

	TCanvas *c10 = new TCanvas("c10","c10",800,800);
	c10->Divide(2,1,0.01,0.01,0);

	c10->cd(1);
	hFingerPMTKE->Draw();
	c10->cd(2);
	hAnaBarPMTKE->Draw();

	c10->DrawClone();
	c10->Print("plots/c10.pdf");

	return c10;

}







