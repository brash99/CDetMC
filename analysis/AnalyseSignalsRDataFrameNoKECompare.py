#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT as root
import numpy as np
import random
from timer import Timer

t = Timer()
t.start()


# In[2]:


print("Environment:")
print("-----------------")

sys_values = True
if (sys_values):
    import sys
    print("sys.executable: ", sys.executable)
    print("-----------------")
    from pprint import pprint as p
    p(sys.path)

print("-----------------")


# In[3]:


fileName1 = "data/AnaBarMC_sieve_v10.root"
fileName2 = "data/AnaBarMC_sieve_v11.root"
treeName = "T"

f1 = root.TFile(fileName1)
f2 = root.TFile(fileName2)
myTree1 = f1.Get(treeName)
myTree2 = f2.Get(treeName)

myTree1.Print()
myTree2.Print()

root.EnableImplicitMT()
d1 = root.RDataFrame(treeName,fileName1)
d2 = root.RDataFrame(treeName,fileName2)

myGeometryData1 = myTree1.GetUserInfo().FindObject("myGeometryData")
myGeometryData2 = myTree2.GetUserInfo().FindObject("myGeometryData")


# In[4]:


#myGeometryData1.Print()
print(myGeometryData1[0])
print(myGeometryData1[0][0])
print(myGeometryData1[0][1])
print(myGeometryData1[0][2])
print(myGeometryData1[0][3])
print(myGeometryData1[1][0])
print(myGeometryData1[1][1])
print(myGeometryData1[1][2])
print(myGeometryData1[1][3])
print(len(myGeometryData1))
type(myGeometryData1)

#myGeometryData2.Print()
print(myGeometryData2[0])
print(myGeometryData2[0][0])
print(myGeometryData2[0][1])
print(myGeometryData2[0][2])
print(myGeometryData2[0][3])
print(myGeometryData2[1][0])
print(myGeometryData2[1][1])
print(myGeometryData2[1][2])
print(myGeometryData2[1][3])
print(len(myGeometryData2))
type(myGeometryData2)




triggerCode = '''
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;
float Photon_min_cut = 80.0;

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

auto fileName = "data/AnaBarMC_9991.root";
auto treeName = "T";

TFile* f = new TFile((TString)fileName,"READ");
TTree* t = (TTree*)f->Get(treeName);

TList* myGeometryData = (TList*)t->GetUserInfo()->FindObject("myGeometryData");

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
        if (ahit) {
            fhit = true;
            trigger = true;
        }
    }
    
    return trigger;
}

int getLayer(int fID) {
    return fID%NUMPADDLE;
}

int getBar(int fID) {
    int iLayer = fID%NUMPADDLE;
    return ((fID-iLayer)/NUMPADDLE)%NUMBARS;
}

int getSide(int fID) {
    int iLayer = fID%NUMPADDLE;
    int iBar = ((fID-iLayer)/NUMPADDLE)%NUMBARS;
    return (((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)%NUMSIDES;
}

int getModule(int fID) {
    int iLayer = fID%NUMPADDLE;
    int iBar = ((fID-iLayer)/NUMPADDLE)%NUMBARS;
    int iSide = (((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)%NUMSIDES;
    return ((((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)/NUMSIDES)%NUMMODULES;
}

int getPlane(int fID) {
    int iLayer = fID%NUMPADDLE;
    int iBar = ((fID-iLayer)/NUMPADDLE)%NUMBARS;
    int iSide = (((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)%NUMSIDES;
    int iModule =  ((((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)/NUMSIDES)%NUMMODULES;
    return (((((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)/NUMSIDES)/NUMMODULES)%NUMLAYERS;
}

float getXOffsetFromTime(int fID, float time) {

    float xoffset;
    int iSide = getSide(fID);
    int iPlane = getPlane(fID);
    TRandom3* fRand = new TRandom3(0);

    /* Old calibration - average time
    if (iPlane == 0) {
        double a = -0.00094362480519633;
        double b = -0.03673192847114299;
        double c = 8.609456854512016;
        double xmin = -25.0;
        double ymin = a*xmin*xmin+b*xmin+c-0.003;
        if (time<ymin) {
            xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
        } else {
            xoffset = -25.0+fRand->Uniform(0.0,12.0);
        }
    } else {
        double a = -0.0009339487175907869;
        double b = -0.03579463613239478;
        double c = 9.59488826868313;
        double xmin = -25.0;
        double ymin = a*xmin*xmin+b*xmin+c-0.003;
        if (time<ymin) {
            xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
        } else {
            xoffset = -25.0+fRand->Uniform(0.0,12.0);
        }
    }
    */
    
    if (iPlane == 0) {
        double a = -0.00037499653844674384;
        double b = -0.051184438607694095;
        double c = 3.5929880196450843;
        double discriminant = b*b - 4*a*(c-time);
        if (discriminant >= 0.0) {
            xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
        } else {
            xoffset = -25.0;
        }
    } else {
        double a = -0.0004090552401677775;
        double b = -0.050453362706664166;
        double c = 4.537007798185197;
        double discriminant = b*b - 4*a*(c-time);
        if (discriminant >= 0.0) {
            xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
        } else {
            xoffset = -25.0;
        }
    }

    if (iSide == 1) {
        xoffset = -xoffset;
    }

    return xoffset;

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

std::vector<double> getDirection(std::vector<double> tp, std::vector<double> bp) {

    std::vector<double> dir;

    double dist = sqrt((tp[0]-bp[0])*(tp[0]-bp[0]) + 
                        (tp[1]-bp[1])*(tp[1]-bp[1]) + 
                         (tp[2]-bp[2])*(tp[2]-bp[2]));
    
    dir.push_back((tp[0]-bp[0])/dist);
    dir.push_back((tp[1]-bp[1])/dist);
    dir.push_back((tp[2]-bp[2])/dist);
    
    return dir;
}

std::vector<double> getPosition(std::vector<std::vector<double>>& sp) {
    //std::cout << "getPosition: " << std::endl;
   
    double x=0;
    double y=0;
    double z=0;
    for (int i = 0; i<sp.size(); i++) {
        std::vector<double> pos = sp[i];
        x += pos[0];
        y += pos[1];
        z += pos[2];
    }
    x = x/sp.size();
    y = y/sp.size();
    z = z/sp.size();
    std::vector<double> tp;
    tp.push_back(x);
    tp.push_back(y);
    tp.push_back(z);
    //cout << tp[0] << " " << tp[1] << " " << tp[2] << std::endl;

    return tp;
}

std::vector<float> getAnaBarPZPMT(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    std::vector<float> v;
    
    std::vector<std::vector<double>> spacePointsTop;
    std::vector<std::vector<double>> spacePointsBottom;

    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                float xdpos, ydpos, zdpos;
                std::vector<double> hitPoint;

                TVectorD* y = (TVectorD*)myGeometryData->At(icount);
                xdpos = (*y)[1]/10.0;
                ydpos = (*y)[2]/10.0;
                zdpos = (*y)[3]/10.0;

                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                //std::cout << "iLayer = " << getLayer(icount) << std::endl;
                //std::cout << "iBar = " << getBar(icount) << std::endl;
                //std::cout << "iSide = " << getSide(icount) << std::endl;
                //std::cout << "iModule = " << getModule(icount) << std::endl;
                //std::cout << "iPlane = " << getPlane(icount) << std::endl;
                pmttime[icount] = PMT_Time[icount];

                float xoffset = getXOffsetFromTime(icount,pmttime[icount]);
                //std::cout << "X offset = " << xoffset << std::endl;
                xdpos = xdpos + xoffset;
                hitPoint.push_back(xdpos);
                hitPoint.push_back(ydpos);
                hitPoint.push_back(zdpos);
                //std::cout << "detector positions: " <<  xdpos << " " << ydpos << " " << zdpos << std::endl;
                
                int iPlane = getPlane(icount);
                if (iPlane == 0) {
                    spacePointsTop.push_back(hitPoint);
                } else {
                    spacePointsBottom.push_back(hitPoint);
                }
            }

        }

        // Fit space points
        std::vector<double> topPosition = getPosition(spacePointsTop);
        std::vector<double> bottomPosition = getPosition(spacePointsBottom);

        std::vector<double> direction = getDirection(topPosition,bottomPosition);

        //std::cout << direction[0] << " " << direction[1] << " " << direction[2] << std::endl;
        
        v.push_back(-direction[2]);

    }
    return v;
}

std::vector<float> getAnaBarPXPMT(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    std::vector<float> v;
    
    std::vector<std::vector<double>> spacePointsTop;
    std::vector<std::vector<double>> spacePointsBottom;

    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                float xdpos, ydpos, zdpos;
                std::vector<double> hitPoint;

                TVectorD* y = (TVectorD*)myGeometryData->At(icount);
                xdpos = (*y)[1]/10.0;
                ydpos = (*y)[2]/10.0;
                zdpos = (*y)[3]/10.0;

                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                //std::cout << "iLayer = " << getLayer(icount) << std::endl;
                //std::cout << "iBar = " << getBar(icount) << std::endl;
                //std::cout << "iSide = " << getSide(icount) << std::endl;
                //std::cout << "iModule = " << getModule(icount) << std::endl;
                //std::cout << "iPlane = " << getPlane(icount) << std::endl;
                pmttime[icount] = PMT_Time[icount];

                float xoffset = getXOffsetFromTime(icount,pmttime[icount]);
                //std::cout << "X offset = " << xoffset << std::endl;
                xdpos = xdpos + xoffset;
                hitPoint.push_back(xdpos);
                hitPoint.push_back(ydpos);
                hitPoint.push_back(zdpos);
                //std::cout << "detector positions: " <<  xdpos << " " << ydpos << " " << zdpos << std::endl;

                int iPlane = getPlane(icount);
                if (iPlane == 0) {
                    spacePointsTop.push_back(hitPoint);
                } else {
                    spacePointsBottom.push_back(hitPoint);
                }
            }

        }

        // Fit space points
        std::vector<double> topPosition = getPosition(spacePointsTop);
        std::vector<double> bottomPosition = getPosition(spacePointsBottom);

        std::vector<double> direction = getDirection(topPosition,bottomPosition);

        //std::cout << direction[0] << " " << direction[1] << " " << direction[2] << std::endl;
        
        v.push_back(-direction[0]);

    }
    return v;
}

std::vector<float> getAnaBarXPMT(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    std::vector<float> v;
    
    std::vector<std::vector<double>> spacePointsTop;
    std::vector<std::vector<double>> spacePointsBottom;

    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                float xdpos, ydpos, zdpos;
                std::vector<double> hitPoint;

                TVectorD* y = (TVectorD*)myGeometryData->At(icount);
                xdpos = (*y)[1]/10.0;
                ydpos = (*y)[2]/10.0;
                zdpos = (*y)[3]/10.0;

                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                //std::cout << "iLayer = " << getLayer(icount) << std::endl;
                //std::cout << "iBar = " << getBar(icount) << std::endl;
                //std::cout << "iSide = " << getSide(icount) << std::endl;
                //std::cout << "iModule = " << getModule(icount) << std::endl;
                //std::cout << "iPlane = " << getPlane(icount) << std::endl;
                pmttime[icount] = PMT_Time[icount];

                float xoffset = getXOffsetFromTime(icount,pmttime[icount]);
                //std::cout << "X offset = " << xoffset << std::endl;
                xdpos = xdpos + xoffset;
                hitPoint.push_back(xdpos);
                hitPoint.push_back(ydpos);
                hitPoint.push_back(zdpos);
                //std::cout << "detector positions: " <<  xdpos << " " << ydpos << " " << zdpos << std::endl;

                int iPlane = getPlane(icount);
                if (iPlane == 0) {
                    spacePointsTop.push_back(hitPoint);
                } else {
                    spacePointsBottom.push_back(hitPoint);
                }
            }

        }

        // Fit space points
        std::vector<double> topPosition = getPosition(spacePointsTop);
        
        v.push_back(topPosition[0]);

    }
    return v;
}

std::vector<float> getAnaBarZPMT(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    std::vector<float> v;
    
    std::vector<std::vector<double>> spacePointsTop;
    std::vector<std::vector<double>> spacePointsBottom;

    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                float xdpos, ydpos, zdpos;
                std::vector<double> hitPoint;

                TVectorD* y = (TVectorD*)myGeometryData->At(icount);
                xdpos = (*y)[1]/10.0;
                ydpos = (*y)[2]/10.0;
                zdpos = (*y)[3]/10.0;

                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                //std::cout << "iLayer = " << getLayer(icount) << std::endl;
                //std::cout << "iBar = " << getBar(icount) << std::endl;
                //std::cout << "iSide = " << getSide(icount) << std::endl;
                //std::cout << "iModule = " << getModule(icount) << std::endl;
                //std::cout << "iPlane = " << getPlane(icount) << std::endl;
                pmttime[icount] = PMT_Time[icount];

                float xoffset = getXOffsetFromTime(icount,pmttime[icount]);
                //std::cout << "X offset = " << xoffset << std::endl;
                xdpos = xdpos + xoffset;
                hitPoint.push_back(xdpos);
                hitPoint.push_back(ydpos);
                hitPoint.push_back(zdpos);
                //std::cout << "detector positions: " <<  xdpos << " " << ydpos << " " << zdpos << std::endl;

                int iPlane = getPlane(icount);
                if (iPlane == 0) {
                    spacePointsTop.push_back(hitPoint);
                } else {
                    spacePointsBottom.push_back(hitPoint);
                }
            }

        }

        // Fit space points
        std::vector<double> topPosition = getPosition(spacePointsTop);
        
        v.push_back(topPosition[2]);

    }
    return v;
}

std::vector<float> getAnaBarPMTTime(bool trigger, int* PMT_Nphotons, float* PMT_Time) {
    
    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                pmttime[icount] = PMT_Time[icount];
                v.push_back(pmttime[icount]);
            }
        }
    }
    return v;
}

std::vector<float> getAnaBarPMTTimeTop(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    std::vector<float> v;

    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                pmttime[icount] = PMT_Time[icount];
                v.push_back(pmttime[icount]);
            }
        }
    }
    return v;
}

std::vector<float> getAnaBarPMTTimeBottom(bool trigger, int* PMT_Nphotons, float* PMT_Time) {

    std::vector<float> v;
   
    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                pmttime[icount] = PMT_Time[icount];
                v.push_back(pmttime[icount]);
            }
        }
    }
    return v;
}


std::vector<float> getAnaBarPMTNPhotons(bool trigger, int* PMT_Nphotons) {
    
    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttot[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];

    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                //std::cout << "getAnaBarPMTNphotons: " << icount << " " << PMT_Nphotons[icount] << std::endl;
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
        for (Int_t icount = Detector_PMT_Offset;icount<Detector_PMT_Offset+4;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
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
            if (PMT_Nphotons[icount]>Photon_min_cut) {
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
            if (PMT_Nphotons[icount]>Photon_min_cut) {
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
'''

root.gInterpreter.Declare(triggerCode)


# In[5]:


fdf1 = d1.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])") \
       .Define("fMass", "getMass(Prim_pdg)") \
       .Define("fMomentum","getMomentum(Prim_E,fMass)") \
       .Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fNewTheta", "getNewTheta(fMomentum,fPy)") \
       .Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)") \
       .Define("trigger2", "getTrigger2(trigger,fNewTheta)") \
       .Define("trigger3", "getTrigger3(trigger,fNewTheta)") \
       .Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)") \
       .Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)") \
       .Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)") \
       .Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)") \
       .Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)") \
       .Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)") \
       .Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)") \
       .Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)") \
       .Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
       .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
       .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])") \
       .Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])") \
       .Define("anaBarXPMT","getAnaBarXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarZPMT","getAnaBarZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPXPMT","getAnaBarPXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPZPMT","getAnaBarPZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTime","getAnaBarPMTTime(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTimeTop","getAnaBarPMTTimeTop(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTimeBottom","getAnaBarPMTTimeBottom(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])") \
       .Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])") \
       .Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])") \
       .Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])") \
       .Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")

triggers = fdf1.Filter("trigger==true").Count()
print('{} entries passed trigger'.format(triggers.GetValue()))

fdft1 = fdf1.Filter("trigger==true")


# In[6]:


fdf2 = d2.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])") \
       .Define("fMass", "getMass(Prim_pdg)") \
       .Define("fMomentum","getMomentum(Prim_E,fMass)") \
       .Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)") \
       .Define("fNewTheta", "getNewTheta(fMomentum,fPy)") \
       .Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)") \
       .Define("trigger2", "getTrigger2(trigger,fNewTheta)") \
       .Define("trigger3", "getTrigger3(trigger,fNewTheta)") \
       .Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)") \
       .Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)") \
       .Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)") \
       .Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)") \
       .Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)") \
       .Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)") \
       .Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)") \
       .Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)") \
       .Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])") \
       .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
       .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
       .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])") \
       .Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])") \
       .Define("anaBarXPMT","getAnaBarXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarZPMT","getAnaBarZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPXPMT","getAnaBarPXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPZPMT","getAnaBarPZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTime","getAnaBarPMTTime(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTimeTop","getAnaBarPMTTimeTop(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarPMTTimeBottom","getAnaBarPMTTimeBottom(trigger,&PMT_Nphotons[0],&PMT_Time[0])") \
       .Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])") \
       .Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])") \
       .Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])") \
       .Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])") \
       .Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")

triggers = fdf2.Filter("trigger==true").Count()
print('{} entries passed trigger'.format(triggers.GetValue()))

fdft2 = fdf2.Filter("trigger==true")


# In[7]:


# Canvas 1

hFingerX1 = fdft1.Histo1D(("h1", "X1", 100, -500.0, 500.0),'fingerXVec')
hFingerY1 = fdft1.Histo1D(("h2", "Y1", 100, -1000.0, 500.0),'fingerYVec')
hFingerZ1 = fdft1.Histo1D(("h3", "Z1", 100, -2000.0, 2000.0),'fingerZVec')
hFingerT1 = fdft1.Histo1D(("h4", "T1", 100, -5.0, 15.0),'fingerTVec')
hFingerX2 = fdft2.Histo1D(("h5", "X2", 100, -500.0, 500.0),'fingerXVec')
hFingerY2 = fdft2.Histo1D(("h6", "Y2", 100, -1000.0, 500.0),'fingerYVec')
hFingerZ2 = fdft2.Histo1D(("h7", "Z2", 100, -2000.0, 2000.0),'fingerZVec')
hFingerT2 = fdft2.Histo1D(("h8", "T2", 100, -5.0, 15.0),'fingerTVec')

c1 = root.TCanvas("c1","c1",800,800)
c1.Divide(2,2,0.01,0.01,0)

c1.cd(1)
hFingerX1.Draw()
hFingerX2.Draw("SAME")
c1.cd(2)
hFingerY1.Draw()
hFingerY2.Draw("SAME")
c1.cd(3)
hFingerZ1.Draw()
hFingerZ2.Draw("SAME")
c1.cd(4)
hFingerT1.Draw()
hFingerT2.Draw("SAME")

print(f"Number of entries GEANT4 v10: {hFingerX1.GetEntries()}")
print(f"Number of entries GEANT4 v11: {hFingerX2.GetEntries()}")

c1.Draw()
c1.Print("plots/c1RA.pdf");


# In[8]:


# Canvas 2

hPrimE1 = fdft1.Histo1D(("h1", "E", 100, 0.0, 7000.0),'Prim_E')
hPrimTh1 = fdft1.Histo1D(("h1", "TH", 100, -0.1, 3.2),'fNewTheta')
hPrimPh1 = fdft1.Histo1D(("h1", "PH", 100, -0.1, 6.4),'fNewPhi')
hPrimPdg1 = fdft1.Histo1D(("h1", "PDG", 100, 0.0, 20.0),'Prim_pdg')
hPrimE2 = fdft2.Histo1D(("h1", "E", 100, 0.0, 7000.0),'Prim_E')
hPrimTh2 = fdft2.Histo1D(("h1", "TH", 100, -0.1, 3.2),'fNewTheta')
hPrimPh2 = fdft2.Histo1D(("h1", "PH", 100, -0.1, 6.4),'fNewPhi')
hPrimPdg2 = fdft2.Histo1D(("h1", "PDG", 100, 0.0, 20.0),'Prim_pdg')

c2 = root.TCanvas("c2","c2",800,800)
c2.Divide(2,2,0.01,0.01,0)

c2.cd(1)
hPrimE2.Draw()
hPrimE1.Draw("SAME")
c2.cd(2)
hPrimTh2.Draw()
hPrimTh1.Draw("SAME")
c2.cd(3)
hPrimPh2.Draw()
hPrimPh1.Draw("SAME")
c2.cd(4)
hPrimPdg2.Draw()
hPrimPdg1.Draw("SAME")

print(f"Number of entries GEANT4 v10: {hPrimE1.GetEntries()}")
print(f"Number of entries GEANT4 v11: {hPrimE2.GetEntries()}")

c2.Draw()
c2.Print("plots/c2RA.pdf")


# In[9]:


def plotSinglePoints(hist):

    x1 = 60.0;
    y1 = 115.0;
    x2 = 70.0;
    y2 = 125.0;
    opacity = 0.9;
    rect7 = root.TBox(x1, y1, x2, y2);
    rect7.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect7);
    x1 = 60.0;
    y1 = -115.0;
    x2 = 70.0;
    y2 = -125.0;
    opacity = 0.9;
    rect8 = root.TBox(x1, y1, x2, y2);
    rect8.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect8);
    x1 = -35.0;
    y1 = -5.0;
    x2 = -45.0;
    y2 = 5.0;
    opacity = 0.9;
    rect9 = root.TBox(x1, y1, x2, y2);
    rect9.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect9);
    x1 = -5.0;
    y1 = 55.0;
    x2 = 5.0;
    y2 = 65.0;
    opacity = 0.9;
    rect10 = root.TBox(x1, y1, x2, y2);
    rect10.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect10);
    x1 = -5.0;
    y1 = -55.0;
    x2 = 5.0;
    y2 = -65.0;
    opacity = 0.9;
    rect11 = root.TBox(x1, y1, x2, y2);
    rect11.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect11);
    x1 = 25.0;
    y1 = -85.0;
    x2 = 35.0;
    y2 = -95.0;
    opacity = 0.9;
    rect12 = root.TBox(x1, y1, x2, y2);
    rect12.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect12);
    x1 = 25.0;
    y1 = 85.0;
    x2 = 35.0;
    y2 = 95.0;
    opacity = 0.9;
    rect13 = root.TBox(x1, y1, x2, y2);
    rect13.SetFillColorAlpha(3, opacity);
    hist.GetListOfFunctions().Add(rect13);

def plotDetector(hist):
    
    opacity=0.2
    x1 = 55.0
    y1 = -61.6
    x2 = -45.2
    y2 = -8.74
    rect1 = root.TBox(x1, y1, x2, y2)
    rect1.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect1)
    x1 = 55.0
    y1 = -8.74
    x2 = -45.2
    y2 = 44.14
    rect2 = root.TBox(x1, y1, x2, y2)
    rect2.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect2)
    x1 = 62.5
    y1 = -114.54
    x2 = -37.7
    y2 = -61.62
    rect3 = root.TBox(x1, y1, x2, y2)
    rect3.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect3)
    x1 = 70.0
    y1 = -167.46
    x2 = -30.2
    y2 = -114.54
    rect4 = root.TBox(x1, y1, x2, y2);
    rect4.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect4);
    x1 = 62.5
    y1 = 44.14
    x2 = -37.7
    y2 = 97.50
    rect5 = root.TBox(x1, y1, x2, y2)
    rect5.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect5)
    x1 = 70.0
    y1 = 97.50
    x2 = -30.2
    y2 = 150.0
    rect6 = root.TBox(x1, y1, x2, y2)
    rect6.SetFillColorAlpha(2, opacity)
    hist.GetListOfFunctions().Add(rect6)

# Canvas 3

hDetectorNhits1 = fdft1.Histo1D(("h1", "Nhits", 50, 0.0, 500.0),'Detector_Nhits')
hDetectorPdg1 = fdft1.Histo1D('anaBarPDG')
hDetectorID1 = fdft1.Histo1D('anaBarID')
hFingerPdg1 = fdft1.Histo1D('fingerPDG')
hFingerID1 = fdft1.Histo1D('fingerID')
hPMTID1 = fdft1.Histo1D('PMT_id')
hAnaBarPMTID1 = fdft1.Histo1D('anaBarPMTID')
hFingerPMTID1 = fdft1.Histo1D('fingerPMTID')

hDetectorNhits2 = fdft2.Histo1D(("h1", "Nhits", 50, 0.0, 500.0),'Detector_Nhits')
hDetectorPdg2 = fdft2.Histo1D('anaBarPDG')
hDetectorID2 = fdft2.Histo1D('anaBarID')
hFingerPdg2 = fdft2.Histo1D('fingerPDG')
hFingerID2 = fdft2.Histo1D('fingerID')
hPMTID2 = fdft2.Histo1D('PMT_id')
hAnaBarPMTID2 = fdft2.Histo1D('anaBarPMTID')
hFingerPMTID2 = fdft2.Histo1D('fingerPMTID')

c3 = root.TCanvas("c3","c3",800,800)
c3.Divide(3,3,0.01,0.01,0)

c3.cd(1)
hDetectorNhits1.Draw()
hDetectorNhits2.Draw("SAME")
c3.cd(2)
hFingerPdg1.Draw()
hFingerPdg2.Draw("SAME")
c3.cd(3)
hDetectorPdg1.Draw()
hDetectorPdg2.Draw("SAME")
c3.cd(4)
hFingerID1.Draw()
hFingerID2.Draw("SAME")
c3.cd(5)
hDetectorID1.Draw()
hDetectorID2.Draw("SAME")
c3.cd(6)
hPMTID1.Draw()
hPMTID2.Draw("SAME")
c3.cd(7)
hFingerPMTID1.Draw()
hFingerPMTID2.Draw("SAME")
c3.cd(8)
hAnaBarPMTID1.Draw()
hAnaBarPMTID2.Draw("SAME")

print(f"Number of entries Detector_Nhits GEANT4 v10: {hDetectorNhits1.GetEntries()}")
print(f"Number of entries Detector_Nhits GEANT4 v11: {hDetectorNhits2.GetEntries()}")
print(f"Number of entries PMT_ID GEANT4 v10: {hPMTID1.GetEntries()}")
print(f"Number of entries PMT_ID GEANT4 v11: {hPMTID2.GetEntries()}")
c3.Draw()
c3.Print("plots/c3RA.pdf")


# In[10]:


# Canvas 33
hAnaBarPMTTime_vs_ID1 = fdft1.Histo2D(("h1", "AnaBar Time vs ID", 100, 0.0, 2500.0,100,0.0,20.0),'anaBarPMTID','anaBarPMTTime')
hAnaBarPMTTime1 = fdft1.Histo1D(("h1", "AnaBar Time", 100, 1.0, 7.0),'anaBarPMTTime')
hAnaBarPMTTimeTop1 = fdft1.Histo1D(("h1", "AnaBar Time Top", 100, 0.0, 13.0),'anaBarPMTTimeTop')
hAnaBarPMTTimeBottom1 = fdft1.Histo1D(("h1", "AnaBar Time Bottom", 100, 0.0, 13.0),'anaBarPMTTimeBottom')
hAnaBarPMTTime_vs_Nphoton1 = fdft1.Histo2D(("h1", "AnaBar Time vs Nphotons", 100, 0.0, 300.0,100,0.0,20.0),'anaBarPMTNPhotons','anaBarPMTTime')
hAnaBarXZPMT1 = fdft1.Histo2D(("h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200),"anaBarXPMT","anaBarZPMT");

hAnaBarPMTTime_vs_ID2 = fdft2.Histo2D(("h1", "AnaBar Time vs ID", 100, 0.0, 2500.0,100,0.0,20.0),'anaBarPMTID','anaBarPMTTime')
hAnaBarPMTTime2 = fdft2.Histo1D(("h1", "AnaBar Time", 100, 1.0, 7.0),'anaBarPMTTime')
hAnaBarPMTTimeTop2 = fdft2.Histo1D(("h1", "AnaBar Time Top", 100, 0.0, 13.0),'anaBarPMTTimeTop')
hAnaBarPMTTimeBottom2 = fdft2.Histo1D(("h1", "AnaBar Time Bottom", 100, 0.0, 13.0),'anaBarPMTTimeBottom')
hAnaBarPMTTime_vs_Nphoton2 = fdft2.Histo2D(("h1", "AnaBar Time vs Nphotons", 100, 0.0, 300.0,100,0.0,20.0),'anaBarPMTNPhotons','anaBarPMTTime')
hAnaBarXZPMT2 = fdft2.Histo2D(("h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200),"anaBarXPMT","anaBarZPMT");

c33 = root.TCanvas("c33","c33",800,1200)
c33.Divide(3,4,0.01,0.01,0)

c33.cd(1)
hAnaBarPMTTime_vs_ID1.Draw("COLZ")
c33.cd(2)
hAnaBarPMTTime_vs_Nphoton1.Draw("COLZ")
c33.cd(3)
hAnaBarPMTTime1.Draw()

c33.cd(4)
hAnaBarPMTTime_vs_ID2.Draw("COLZ")
c33.cd(5)
hAnaBarPMTTime_vs_Nphoton2.Draw("COLZ")
c33.cd(6)
hAnaBarPMTTime2.Draw()

c33.cd(7)
hAnaBarPMTTimeTop1.Draw()
c33.cd(8)
hAnaBarPMTTimeBottom1.Draw()
c33.cd(9)
hAnaBarXZPMT1.Draw("COLZ")
plotDetector(hAnaBarXZPMT1)
plotSinglePoints(hAnaBarXZPMT1)

c33.cd(10)
hAnaBarPMTTimeTop2.Draw()
c33.cd(11)
hAnaBarPMTTimeBottom2.Draw()
c33.cd(12)
hAnaBarXZPMT2.Draw("COLZ")
plotDetector(hAnaBarXZPMT2)
plotSinglePoints(hAnaBarXZPMT2)

print(f"Number of entries PMT Time GEANT4 v10: {hAnaBarPMTTime1.GetEntries()}")
print(f"Number of entries PMT Time GEANT4 v11: {hAnaBarPMTTime2.GetEntries()}")

c33.Draw()
c33.Print("plots/c33.pdf")


# In[12]:


hAnaBarXZPMT1 = fdft1.Histo2D(("h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200),"anaBarXPMT","anaBarZPMT")
hAnaBarXZPMT2 = fdft2.Histo2D(("h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200),"anaBarXPMT","anaBarZPMT")

c34 = root.TCanvas("c34","c34",800,1600)
c34.Divide(2,2,0.01,0.01,0)

c34.cd(1)
hAnaBarXZPMT1.Draw("SURF1")

c34.cd(2)
hAnaBarXZPMT1.Draw("CONT1")
plotDetector(hAnaBarXZPMT1)
#plotSinglePoints(hAnaBarXZPMT);

c34.cd(3)
hAnaBarXZPMT2.Draw("SURF1")

c34.cd(4)
hAnaBarXZPMT2.Draw("CONT1")
plotDetector(hAnaBarXZPMT2)
#plotSinglePoints(hAnaBarXZPMT);
                  
c34.Print("plots/c34.pdf");


# In[13]:


# Canvas 4

hFingerEd1 = fdft1.Histo1D(("h1", "Finger EDep", 100, 0.0, 10.0),'fingerEd')
hFingerPMTNphot1 = fdft1.Histo1D(("h1", "Finger Npe", 100, 0.0, 20.0),'fingerPMTNPhotons')
hAnaBarPMTNphot1 = fdft1.Histo1D(("h1", "AnaBar Npe", 100, 0.0, 200.0),'anaBarPMTNPhotons')
hAnaBarEd1 = fdft1.Histo1D(("h1", "AnaBar Edep", 100, 0.0, 10.0),'anaBarEd')

hFingerEd2 = fdft2.Histo1D(("h1", "Finger EDep", 100, 0.0, 10.0),'fingerEd')
hFingerPMTNphot2 = fdft2.Histo1D(("h1", "Finger Npe", 100, 0.0, 20.0),'fingerPMTNPhotons')
hAnaBarPMTNphot2 = fdft2.Histo1D(("h1", "AnaBar Npe", 100, 0.0, 200.0),'anaBarPMTNPhotons')
hAnaBarEd2 = fdft2.Histo1D(("h1", "AnaBar Edep", 100, 0.0, 10.0),'anaBarEd')

c4 = root.TCanvas("c4","c4",800,800)

c4.cd()
pad = root.TPad("pad","pad",0.01,0.51,0.50,0.99)
pad.Draw()
pad.cd()
hFingerEd2.GetXaxis().SetRangeUser(1.0,10)
hFingerEd2.Draw();
hFingerEd1.GetXaxis().SetRangeUser(1.0,10)
hFingerEd1.Draw("SAME");

c4.cd()
pad = root.TPad("pad","pad",0.51,0.51,0.99,0.99)
pad.Draw()
pad.cd()
hFingerPMTNphot2.GetXaxis().SetRangeUser(-10,250)
hFingerPMTNphot2.Draw()
hFingerPMTNphot1.GetXaxis().SetRangeUser(-10,250)
hFingerPMTNphot1.Draw("SAME")

c4.cd()
pad = root.TPad("pad","pad",0.01,0.01,0.50,0.50)
#pad.SetLogy()
pad.Draw()
pad.cd()
hAnaBarEd2.GetXaxis().SetRangeUser(1.0,10)
hAnaBarEd2.Draw();
hAnaBarEd1.GetXaxis().SetRangeUser(1.0,10)
hAnaBarEd1.Draw("SAME");

c4.cd()
pad = root.TPad("pad","pad",0.51,0.01,0.99,0.50)
#pad.SetLogy()
pad.Draw()
pad.cd()
hAnaBarPMTNphot2.GetXaxis().SetRangeUser(-20,180)
hAnaBarPMTNphot2.Draw()
hAnaBarPMTNphot1.GetXaxis().SetRangeUser(-20,180)
hAnaBarPMTNphot1.Draw("SAME")

print(f"Number of entries AnaBar PMT Number of Photons GEANT4 v10: {hAnaBarPMTNphot1.GetEntries()}")
print(f"Number of entries AnaBar PMT Number of Photons GEANT4 v11: {hAnaBarPMTNphot2.GetEntries()}")

c4.Draw()
c4.Print("plots/c4RA.pdf")


# In[14]:


# Canvas 5

hAnaBarX1 = fdft1.Histo1D(("h1", "AnaBar X", 1000, -800.0, 800.0),'anaBarXVec')
hAnaBarY1 = fdft1.Histo1D(("h1", "AnaBar Y", 1000, -400.0, 100.0),'anaBarYVec')
hAnaBarZ1 = fdft1.Histo1D(("h1", "AnaBar Z", 1000, -2000.0, 2000.0),'anaBarZVec')
hAnaBarT1 = fdft1.Histo1D(("h1", "AnaBar T", 1000, 0.0, 10.0),'anaBarTVec')
hAnaBarX2 = fdft2.Histo1D(("h1", "AnaBar X", 1000, -800.0, 800.0),'anaBarXVec')
hAnaBarY2 = fdft2.Histo1D(("h1", "AnaBar Y", 1000, -400.0, 100.0),'anaBarYVec')
hAnaBarZ2 = fdft2.Histo1D(("h1", "AnaBar Z", 1000, -2000.0, 2000.0),'anaBarZVec')
hAnaBarT2 = fdft2.Histo1D(("h1", "AnaBar T", 1000, 0.0, 10.0),'anaBarTVec')

c5 = root.TCanvas("c5","c5",800,800)
c5.Divide(2,2,0.01,0.01,0)

c5.cd(1)
hAnaBarX2.Draw()
c5.cd(2)
hAnaBarY2.Draw()
c5.cd(3)
hAnaBarZ2.Draw()
c5.cd(4)
hAnaBarT2.Draw()
c5.cd(1)
hAnaBarX1.Draw("SAME")
c5.cd(2)
hAnaBarY1.Draw("SAME")
c5.cd(3)
hAnaBarZ1.Draw("SAME")
c5.cd(4)
hAnaBarT1.Draw("SAME")


print(f"Number of entries AnaBarX1 GEANT4 v10: {hAnaBarX1.GetEntries()}")
print(f"Number of entries AnaBarX2 GEANT4 v11: {hAnaBarX2.GetEntries()}")
print(f"Number of entries AnaBarY1 GEANT4 v10: {hAnaBarY1.GetEntries()}")
print(f"Number of entries AnaBarY2 GEANT4 v11: {hAnaBarY2.GetEntries()}")
print(f"Number of entries AnaBarZ1 GEANT4 v10: {hAnaBarZ1.GetEntries()}")
print(f"Number of entries AnaBarZ2 GEANT4 v11: {hAnaBarZ2.GetEntries()}")
print(f"Number of entries AnaBarT1 GEANT4 v10: {hAnaBarT1.GetEntries()}")
print(f"Number of entries AnaBarT2 GEANT4 v11: {hAnaBarT2.GetEntries()}")

c5.Draw()
c5.Print("plots/c5RA.pdf");


# In[ ]:





# In[15]:


hE1vsE21 = fdft1.Histo2D(("h2", "E1 vs E2", 100, 0.01, 3.0, 100, 0.01, 30.0),"fingerEd","anaBarEdTotal")
hE1vsE22 = fdft2.Histo2D(("h2", "E1 vs E2", 100, 0.01, 3.0, 100, 0.01, 30.0),"fingerEd","anaBarEdTotal")

c6 = root.TCanvas("c6", "c6", 800, 800)
c6.Divide(2,1, 0.01, 0.01, 0)

c6.cd(1)
hE1vsE21.Draw("COLZ")
c6.cd(2)
hE1vsE22.Draw("COLZ")

c6.Draw()
c6.Print("plots/c6RA.pdf")


# In[16]:


hFinger_Edep_vs_Nphot1 = fdft1.Filter("trigger2").Histo2D(("h3", "Finger Edep vs Nphot", 100, 0.01, 200.0, 100, 0.01, 3.0),"fingerPMTNPhotons","fingerEd")
hAnaBar_Edep_vs_Nphot1 = fdft1.Filter("trigger2").Histo2D(("h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0),"anaBarEdTotal","anaBarNPhotonsTotal")
hNphot0_vs_Nphot1 = fdft1.Filter("trigger2").Histo2D(("h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 200.0),"anaBarNPhotonsTotal","fingerPMTNPhotons")
hFinger_Edep_vs_Nphot2 = fdft2.Filter("trigger2").Histo2D(("h3", "Finger Edep vs Nphot", 100, 0.01, 200.0, 100, 0.01, 3.0),"fingerPMTNPhotons","fingerEd")
hAnaBar_Edep_vs_Nphot2 = fdft2.Filter("trigger2").Histo2D(("h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0),"anaBarEdTotal","anaBarNPhotonsTotal")
hNphot0_vs_Nphot2 = fdft2.Filter("trigger2").Histo2D(("h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 200.0),"anaBarNPhotonsTotal","fingerPMTNPhotons")

c7 = root.TCanvas("c7", "c7", 800, 800)
c7.Divide(4,2, 0.01, 0.01, 0)

#c7PE_MeV = root.TCanvas("c7PE_MeV", "c7PE_MeV", 800,800)
#c7Profile = root.TCanvas("c7Profile", "c7Profile", 800,800)

c7.cd(1)
hFinger_Edep_vs_Nphot1.Draw("COLZ")
c7.cd(2)
hAnaBar_Edep_vs_Nphot1.Draw("COLZ")
c7.cd(3)
hNphot0_vs_Nphot1.Draw("COLZ")
c7.cd(4)
prof = hAnaBar_Edep_vs_Nphot1.ProfileX()
prof.Fit("pol1")
c7.cd(5)
hFinger_Edep_vs_Nphot2.Draw("COLZ")
c7.cd(6)
hAnaBar_Edep_vs_Nphot2.Draw("COLZ")
c7.cd(7)
hNphot0_vs_Nphot2.Draw("COLZ")
c7.cd(8)
prof = hAnaBar_Edep_vs_Nphot2.ProfileX()
prof.Fit("pol1")


# In[17]:


c7.Draw()
c7.Print("plots/c7RA.pdf")


# In[18]:


hFinger_Edep_vs_Nphot1b = fdft1.Filter("trigger3").Histo2D(("h3", "Finger Edep vs Nphot", 100, 0.01, 200.0, 100, 0.01, 3.0),"fingerPMTNPhotons","fingerEd")
hAnaBar_Edep_vs_Nphot1b = fdft1.Filter("trigger3").Histo2D(("h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0),"anaBarEdTotal","anaBarNPhotonsTotal")
hNphot0_vs_Nphot1b = fdft1.Filter("trigger3").Histo2D(("h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 200.0),"anaBarNPhotonsTotal","fingerPMTNPhotons")
hFinger_Edep_vs_Nphot2b = fdft2.Filter("trigger3").Histo2D(("h3", "Finger Edep vs Nphot", 100, 0.01, 200.0, 100, 0.01, 3.0),"fingerPMTNPhotons","fingerEd")
hAnaBar_Edep_vs_Nphot2b = fdft2.Filter("trigger3").Histo2D(("h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0),"anaBarEdTotal","anaBarNPhotonsTotal")
hNphot0_vs_Nphot2b = fdft2.Filter("trigger3").Histo2D(("h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 200.0),"anaBarNPhotonsTotal","fingerPMTNPhotons")

c8 = root.TCanvas("c8", "c8", 800, 800)
c8.Divide(4,2, 0.01, 0.01, 0)

#c7PE_MeV = root.TCanvas("c7PE_MeV", "c7PE_MeV", 800,800)
#c7Profile = root.TCanvas("c7Profile", "c7Profile", 800,800)

c8.cd(1)
hFinger_Edep_vs_Nphot1b.Draw("COLZ")
c8.cd(2)
hAnaBar_Edep_vs_Nphot1b.Draw("COLZ")
c8.cd(3)
hNphot0_vs_Nphot1b.Draw("COLZ")
c8.cd(4)
prof = hAnaBar_Edep_vs_Nphot1b.ProfileX()
prof.Fit("pol1")
c8.cd(5)
hFinger_Edep_vs_Nphot2b.Draw("COLZ")
c8.cd(6)
hAnaBar_Edep_vs_Nphot2b.Draw("COLZ")
c8.cd(7)
hNphot0_vs_Nphot2b.Draw("COLZ")
c8.cd(8)
prof = hAnaBar_Edep_vs_Nphot2b.ProfileX()
prof.Fit("pol1")


# In[19]:


c8.Draw()
c8.Print("plots/c8RA.pdf")


# In[20]:


#NUMPADDLE=14
#
#hAnaBarEdAll = []
#
#for i in range(NUMPADDLE):
#    name = ("AnaBarEd%d" % i)
#    title = ("AnaBar Energy Deposited A%d" % i)
#    name2 = ("anaBarEd[%d]" %i)
#    hAnaBarEdAll.append(fdft.Define(name,name2).Filter("trigger3").Histo1D((name, title, 100, 0.01, 10.0),name))
#    
#hAnaBarEdAllCut = []
#
#for i in range(NUMPADDLE):
#    name = ("AnaBarEd%dCut" % i)
#    title = ("AnaBar Energy Deposited A%d" % i)
#    name2 = ("anaBarEd[%d]" %i)
#    name3 = ("trigger3 && anaBarPMTNPhotons[%d]>100.0" % i)
#    hAnaBarEdAllCut.append(fdft.Define(name,name2).Filter(name3).Histo1D((name, title, 100, 0.01, 10.0),name))
#
#cEd = root.TCanvas("cEd", "cEd", 800,800)
#cEd.Divide(4,4)
#
#means = []
#meanErr = []
#
#start = 5.5
#gf = root.TF1("gf", "gaus", start, 10.0)
#
#for i in range(NUMPADDLE):
#    
#    print ("Paddle = ",i+1)
#
#    cEd.cd(i+1)
#    
#    hAnaBarEdAll[i].Draw()
#    hAnaBarEdAllCut[i].SetLineColor(2)
#    hAnaBarEdAllCut[i].Draw("SAME")
#    
#cEd.Draw()
#cEd.Print("plots/cEdRA.pdf")


# In[21]:


#hAnaBarPMTNphotArray = []
#
#for i in range(NUMPADDLE):  
#    name = ("AnaBarPMTNphotA%d" % i)
#    title = ("AnaBar_PMT_Number_of_Photons_A%d" % i)
#    name2 = ("anaBarPMTNPhotons[%d]" % i)
#    hAnaBarPMTNphotArray.append(fdft.Define(name,name2).Filter("trigger").Histo1D((name, title, 200, -20, 180.0),name))
#    
#hAnaBarPMTNphotCut = []
#    
#for i in range(NUMPADDLE):  
#    name = ("AnaBarPMTNphotA%dCut" % i)
#    title = ("AnaBar_PMT_Number_of_Photons_A%d_Cut" % i)
#    name2 = ("anaBarPMTNPhotons[%d]" % i)
#    hAnaBarPMTNphotCut.append(fdft.Define(name,name2).Filter("trigger3").Histo1D((name, title, 200, -20, 180.0),name))
#    
#hAnaBarPMTNoiseCutNphot = []
#
#for i in range(NUMPADDLE):
#    name = ("AnaBarPMTNoiseCutNphotA%d" % i)
#    title = ("AnaBar PMT Number of Photons A%d" % i)
#    name2 = ("anaBarPMTNPhotons[%d]" % i)
#    name3 = ("trigger3 && anaBarEd[%d]>4.0" % i)
#    hAnaBarPMTNoiseCutNphot.append(fdft.Define(name,name2).Filter(name3).Histo1D((name, title, 200, -20, 180.0),name))
#
#c9 = root.TCanvas("c9", "c9", 800,800)
#
#print("Fitting ...\n")
#fr = [float, float]
#fp, fpe = [float, float, float, float], [float, float, float, float]
#pllo = [0.05, 0.5, 1.0, 0.04]
#plhi = [10.0, 50.0, 10000.0, 5.0]
#sv = [1.8, 5.0, 1400.0, 3.0]
#chisqr = float
#ndf = int
#SNRPeak, SNRFWHM = float, float
#
#
#for i in range(NUMPADDLE):
#        c9.cd()
#        xl = 0.25*(i%4)
#        xh = 0.25*(i%4)+0.25
#        yl = 0.75 - 0.25*int(i/4)
#        yh = 0.75 - 0.25*int(i/4) + 0.25
#        #print (i,xl,yl,xh,yh)
#        pad = root.TPad("pad","pad",xl,yl,xh,yh)
#        pad.SetLogy(True)
#        pad.Draw()
#        pad.cd()
#        hAnaBarPMTNphotArray[i].Draw()
#        hAnaBarPMTNoiseCutNphot[i].SetLineColor(2)
#        hAnaBarPMTNoiseCutNphot[i].Draw("SAME")
#        #hAnaBarPMTNoiseCutNphot[i].Fit("gaus")


# In[22]:


#c9.Draw()
#c9.Print("plots/c9RA.pdf")


# In[ ]:





# In[23]:


hAnaBarMult1 = fdft1.Histo1D(("h1", "CDet Multiplicity", 20, 0, 20),"imult")
hAnaBarMult2 = fdft2.Histo1D(("h1", "CDet Multiplicity", 20, 0, 20),"imult")

c11 = root.TCanvas("c11", "c11", 800,800)
c11.Divide(1,1, 0.01, 0.01, 0)

c11.cd(1)
hAnaBarMult1.Draw()
hAnaBarMult2.Draw("SAME")

c11.Draw()
c11.Print("plots/c11RA.pdf")


# In[24]:


hPrimPx1 = fdft1.Histo1D(("h1", "Px", 100, -1000, 1000),"fPx")
hPrimPy1 = fdft1.Histo1D(("h1", "Px", 100, -8000, 1000),"fPy")
hPrimPz1 = fdft1.Histo1D(("h1", "Px", 100, -2000, 2000),"fPz")
hPrimPx2 = fdft2.Histo1D(("h1", "Px", 100, -1000, 1000),"fPx")
hPrimPy2 = fdft2.Histo1D(("h1", "Px", 100, -8000, 1000),"fPy")
hPrimPz2 = fdft2.Histo1D(("h1", "Px", 100, -2000, 2000),"fPz")

c12 = root.TCanvas("c12", "c12", 800,800)
c12.Divide(2,2, 0.01, 0.01, 0)

c12.cd(1)
hPrimPx2.Draw()
c12.cd(2)
hPrimPy2.Draw()
c12.cd(3)
hPrimPz2.Draw()
c12.cd(1)
hPrimPx1.Draw("SAME")
c12.cd(2)
hPrimPy1.Draw("SAME")
c12.cd(3)
hPrimPz1.Draw("SAME")

c12.Draw()
c12.Print("plots/c12RA.pdf")


# In[ ]:





# In[25]:


hPx_vs_x1 = fdft1.Filter("trigger2").Histo2D(("h33", "Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0),"anaBarXVec","fPx")
hPz_vs_z1 = fdft1.Filter("trigger2").Histo2D(("h34", "Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0),"anaBarZVec","fPz")
hz_vs_x1 = fdft1.Filter("trigger2").Histo2D(("h35", "z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0),"anaBarXVec","anaBarZVec")
hPx_vs_x2 = fdft2.Filter("trigger2").Histo2D(("h33", "Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0),"anaBarXVec","fPx")
hPz_vs_z2 = fdft2.Filter("trigger2").Histo2D(("h34", "Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0),"anaBarZVec","fPz")
hz_vs_x2 = fdft2.Filter("trigger2").Histo2D(("h35", "z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0),"anaBarXVec","anaBarZVec")

hPrimXZ1 = fdft1.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z")
hPrimXZ2 = fdft2.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z")

c13 = root.TCanvas("c13", "c13", 800, 800)
c13.Divide(4,2, 0.01, 0.01, 0)

c13.cd(1)
hPx_vs_x1.Draw("COLZ")
c13.cd(2)
hPz_vs_z1.Draw("COLZ")
c13.cd(3)
hz_vs_x1.Draw("COLZ")
c13.cd(4)
plotDetector(hPrimXZ1);
hPrimXZ1.Draw('COLZ')

c13.cd(5)
hPx_vs_x2.Draw("COLZ")
c13.cd(6)
hPz_vs_z2.Draw("COLZ")
c13.cd(7)
hz_vs_x2.Draw("COLZ")
c13.cd(8)
plotDetector(hPrimXZ2);
hPrimXZ2.Draw('COLZ')


c13.Draw()
c13.Print("plots/c13RA.pdf")


# In[26]:


hPrimX1 = fdft1.Histo1D(("h99","X_vtx", 100, -80,80),"Prim_X")
hPrimY1 = fdft1.Histo1D(("h99","Y_vtx", 100, 0,40),"Prim_Y")
hPrimZ1 = fdft1.Histo1D(("h99","Z_vtx", 100, -200,200),"Prim_Z")
hPrimXZ1 = fdft1.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z")
hPrimX2 = fdft2.Histo1D(("h99","X_vtx", 100, -80,80),"Prim_X")
hPrimY2 = fdft2.Histo1D(("h99","Y_vtx", 100, 0,40),"Prim_Y")
hPrimZ2 = fdft2.Histo1D(("h99","Z_vtx", 100, -200,200),"Prim_Z")
hPrimXZ2 = fdft2.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z")

c14 = root.TCanvas("c14","c14",800,800)
c14.Divide(3,2,0.01,0.01,0)

c14.cd(1);
hPrimX1.Draw();
c14.cd(2); 
hPrimY1.Draw();
c14.cd(3);
hPrimZ1.Draw();
c14.cd(4);
hPrimXZ1.Draw("COLZ");
c14.cd(1);
hPrimX2.Draw("SAME");
c14.cd(2);
hPrimY2.Draw("SAME");
c14.cd(3);
hPrimZ2.Draw("SAME");
c14.cd(5);
hPrimXZ2.Draw("COLZ");

c14.Draw();
c14.Print("plots/c14.pdf");


# In[27]:


hPx_vs_x1 = fdft1.Filter("trigger2").Histo2D(("h33", "Px vs x", 100, -80.0, 80.0, 100, -1.0, 1.0),"anaBarXPMT","anaBarPXPMT");
hPz_vs_z1 = fdft1.Filter("trigger2").Histo2D(("h34", "Pz vs z", 100, -240.0, 240.0, 100, -1.0, 1.0),"anaBarZPMT","anaBarPZPMT");
hz_vs_x1 = fdft1.Filter("trigger2").Histo2D(("h35", "z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0),"anaBarXPMT","anaBarZPMT");
hPrimXZ1 = fdft1.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z");
hPx_vs_x2 = fdft2.Filter("trigger2").Histo2D(("h33", "Px vs x", 100, -80.0, 80.0, 100, -1.0, 1.0),"anaBarXPMT","anaBarPXPMT");
hPz_vs_z2 = fdft2.Filter("trigger2").Histo2D(("h34", "Pz vs z", 100, -240.0, 240.0, 100, -1.0, 1.0),"anaBarZPMT","anaBarPZPMT");
hz_vs_x2 = fdft2.Filter("trigger2").Histo2D(("h35", "z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0),"anaBarXPMT","anaBarZPMT");
hPrimXZ2 = fdft2.Histo2D(("h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0),"Prim_X","Prim_Z");



c15 = root.TCanvas("c15","c15",800,800);
c15.Divide(4,2,0.01,0.01,0);

c15.cd(1);
hPx_vs_x1.Draw("COLZ");
c15.cd(2);
hPz_vs_z1.Draw("COLZ");
c15.cd(3);
hz_vs_x1.Draw("COLZ");
plotDetector(hz_vs_x1);
c15.cd(4);
hPrimXZ1.Draw("COLZ");
plotDetector(hPrimXZ1);

c15.cd(5);
hPx_vs_x2.Draw("COLZ");
c15.cd(6);
hPz_vs_z2.Draw("COLZ");
c15.cd(7);
hz_vs_x2.Draw("COLZ");
plotDetector(hz_vs_x2);
c15.cd(8);
hPrimXZ2.Draw("COLZ");
plotDetector(hPrimXZ2);

c15.Draw();
c15.Print("plots/c15.pdf");


# In[28]:


hx_vs_x1 = fdft1.Filter("trigger2").Histo2D(("h33", "X vs XPMT", 100, -80.0, 80.0, 100, -80.0, 80.0),"anaBarXPMT","Prim_X");
hz_vs_z1 = fdft1.Filter("trigger2").Histo2D(("h34", "Z vs ZPMT", 100, -240.0, 240.0, 100, -240.0, 240.0),"anaBarZPMT","Prim_Z");
hPx_vs_Px1 = fdft1.Filter("trigger2").Histo2D(("h35", "Px vs PxPMT", 100, -1.0, 1.0, 100, -800.0, 800.0),"anaBarPXPMT","fPx");
hPz_vs_Pz1 = fdft1.Filter("trigger2").Histo2D(("h99", "Pz vs PzPMT", 100, -0.6, 0.6, 100, -2400.0, 2400.0),"anaBarPZPMT","fPz");
hx_vs_x2 = fdft2.Filter("trigger2").Histo2D(("h33", "X vs XPMT", 100, -80.0, 80.0, 100, -80.0, 80.0),"anaBarXPMT","Prim_X");
hz_vs_z2 = fdft2.Filter("trigger2").Histo2D(("h34", "Z vs ZPMT", 100, -240.0, 240.0, 100, -240.0, 240.0),"anaBarZPMT","Prim_Z");
hPx_vs_Px2 = fdft2.Filter("trigger2").Histo2D(("h35", "Px vs PxPMT", 100, -1.0, 1.0, 100, -800.0, 800.0),"anaBarPXPMT","fPx");
hPz_vs_Pz2 = fdft2.Filter("trigger2").Histo2D(("h99", "Pz vs PzPMT", 100, -0.6, 0.6, 100, -2400.0, 2400.0),"anaBarPZPMT","fPz");

c16 = root.TCanvas("c16","c16",800,800);
c16.Divide(4,2,0.01,0.01,0);

c16.cd(1);
hx_vs_x1.Draw("COLZ");
c16.cd(2);
hz_vs_z1.Draw("COLZ");
c16.cd(3);
hPx_vs_Px1.Draw("COLZ");
c16.cd(4);
hPz_vs_Pz1.Draw("COLZ");
c16.cd(5);
hx_vs_x2.Draw("COLZ");
c16.cd(6);
hz_vs_z2.Draw("COLZ");
c16.cd(7);
hPx_vs_Px2.Draw("COLZ");
c16.cd(8);
hPz_vs_Pz2.Draw("COLZ");

c16.Draw();
c16.Print("plots/c16.pdf");


# In[29]:


t.stop()


# In[ ]:





# In[ ]:




