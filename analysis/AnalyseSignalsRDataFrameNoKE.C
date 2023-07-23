#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>

using namespace std;
using RNode = ROOT::RDF::RNode;

std::vector<RNode> v;
TList* myGeometryData;
int global_run_number;
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;
float Photon_min_cut = 75.0;

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
		xoffset -25.0;
	}
    } else {
        double a = -0.0004090552401677775;
        double b = -0.050453362706664166;
        double c = 4.537007798185197;
	double discriminant = b*b - 4*a*(c-time);
	if (discriminant >= 0.0) {
        	xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
	} else {
		xoffset -25.0;
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
    
    float pmttime[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];

    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
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

//RNode AnalyseSignalsRDataFrameNoKE(int run_number = 4000) {
void AnalyseSignalsRDataFrameNoKE(int run_number = 4000) {

	global_run_number = run_number;
	std::cout << run_number << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);
	auto fileName = "data/AnaBarMC_"+std::to_string(run_number)+".root";
	auto treeName = "T";

        TFile* f = new TFile((TString)fileName,"READ");
        TTree* t = (TTree*)f->Get(treeName);

	ROOT::RDataFrame d(treeName,fileName);

        myGeometryData = (TList*)t->GetUserInfo()->FindObject("myGeometryData");

        //myGeometryData->Print();

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
                        .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarXPMT","getAnaBarXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarZPMT","getAnaBarZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPXPMT","getAnaBarPXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPZPMT","getAnaBarPZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTime","getAnaBarPMTTime(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTimeTop","getAnaBarPMTTimeTop(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTimeBottom","getAnaBarPMTTimeBottom(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");

	//auto entries2 = fdf.Count();
	//cout << *entries2 << " entries in Expanded Dataframe with no filter" << endl;

	auto triggers = fdf.Filter("trigger==true").Count();
	cout << *triggers << " entries passed Main trigger" << endl;

        auto fdft = fdf.Filter("trigger==true");
        v.push_back(fdft);

	//return fdft;
}

TCanvas* plotC1(){

  //RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

  auto hFingerX = v[0].Histo1D("fingerXVec");
  auto hFingerY = v[0].Histo1D("fingerYVec");
  auto hFingerZ = v[0].Histo1D("fingerZVec");
  auto hFingerT = v[0].Histo1D("fingerTVec");
  
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

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimE = v[0].Histo1D("Prim_E");
	auto hPrimTh = v[0].Histo1D("fNewTheta");
	auto hPrimPh = v[0].Histo1D("fNewPhi");
	auto hPrimPdg = v[0].Histo1D("Prim_pdg");

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

void plotDetector(ROOT::RDF::RResultPtr<TH2D> hist) {

    double opacity=0.2;
    double x1 = 55.0;
    double y1 = -61.6;
    double x2 = -45.2;
    double y2 = -8.74;
    TBox *rect1 = new TBox(x1, y1, x2, y2);
    rect1->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect1);
    x1 = 55.0;
    y1 = -8.74;
    x2 = -45.2;
    y2 = 44.14;
    TBox *rect2 = new TBox(x1, y1, x2, y2);
    rect2->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect2);
    x1 = 62.5;
    y1 = -114.54;
    x2 = -37.7;
    y2 = -61.62;
    TBox *rect3 = new TBox(x1, y1, x2, y2);
    rect3->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect3);
    x1 = 70.0;
    y1 = -167.46;
    x2 = -30.2;
    y2 = -114.54;
    TBox *rect4 = new TBox(x1, y1, x2, y2);
    rect4->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect4);
    x1 = 62.5;
    y1 = 44.14;
    x2 = -37.7;
    y2 = 97.50;
    TBox *rect5 = new TBox(x1, y1, x2, y2);
    rect5->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect5);
    x1 = 70.0;
    y1 = 97.50;
    x2 = -30.2;
    y2 = 150.0;
    TBox *rect6 = new TBox(x1, y1, x2, y2);
    rect6->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect6);

}

void plotSinglePoints(ROOT::RDF::RResultPtr<TH2D> hist) {

    double x1 = 60.0;
    double y1 = 115.0;
    double x2 = 70.0;
    double y2 = 125.0;
    double opacity = 0.9;
    TBox *rect7 = new TBox(x1, y1, x2, y2);
    rect7->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect7);
    x1 = 60.0;
    y1 = -115.0;
    x2 = 70.0;
    y2 = -125.0;
    opacity = 0.9;
    TBox *rect8 = new TBox(x1, y1, x2, y2);
    rect8->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect8);
    x1 = -35.0;
    y1 = -5.0;
    x2 = -45.0;
    y2 = 5.0;
    opacity = 0.9;
    TBox *rect9 = new TBox(x1, y1, x2, y2);
    rect9->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect9);
    x1 = -5.0;
    y1 = 55.0;
    x2 = 5.0;
    y2 = 65.0;
    opacity = 0.9;
    TBox *rect10 = new TBox(x1, y1, x2, y2);
    rect10->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect10);
    x1 = -5.0;
    y1 = -55.0;
    x2 = 5.0;
    y2 = -65.0;
    opacity = 0.9;
    TBox *rect11 = new TBox(x1, y1, x2, y2);
    rect11->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect11);
    x1 = 25.0;
    y1 = -85.0;
    x2 = 35.0;
    y2 = -95.0;
    opacity = 0.9;
    TBox *rect12 = new TBox(x1, y1, x2, y2);
    rect12->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect12);
    x1 = 25.0;
    y1 = 85.0;
    x2 = 35.0;
    y2 = 95.0;
    opacity = 0.9;
    TBox *rect13 = new TBox(x1, y1, x2, y2);
    rect13->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect13);

}


TCanvas* plotC33(){

    auto hAnaBarPMTTime_vs_ID = v[0].Histo2D({"h1", "AnaBar Time vs ID", 100, 0.0, 2500.0,100,0.0,20.0},"anaBarPMTID","anaBarPMTTime");
    auto hAnaBarXZPMT = v[0].Histo2D({"h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200},"anaBarXPMT","anaBarZPMT");
    auto hAnaBarPMTTime = v[0].Histo1D({"h1","AnaBar Time",100,0.0,20.0},"anaBarPMTTime");
    auto hAnaBarPMTTimeTop = v[0].Histo1D({"h1","AnaBar Time Top",100,0.0,20.0},"anaBarPMTTimeTop");
    auto hAnaBarPMTTimeBottom = v[0].Histo1D({"h1","AnaBar Time Bottom",100,0.0,20.0},"anaBarPMTTimeBottom");
    auto hAnaBarPMTTime_vs_Nphoton = v[0].Histo2D({"h1", "AnaBar Time vs Nphotons", 100, 0.0, 300.0,100,0.0,20.0},"anaBarPMTNPhotons","anaBarPMTTime");

    TCanvas* c33 = new TCanvas("c33","c33",800,800);
    c33->Divide(3,2,0.01,0.01,0);

    c33->cd(1);
    hAnaBarPMTTime_vs_ID->Draw("COLZ");
    c33->cd(2);
    hAnaBarPMTTime_vs_Nphoton->Draw("COLZ");
    c33->cd(3);
    hAnaBarPMTTime->Draw();
    hAnaBarPMTTimeTop->Draw("SAME");
    hAnaBarPMTTimeBottom->Draw("SAME");
    c33->cd(4);
    hAnaBarPMTTimeTop->Draw();
    c33->cd(5);
    hAnaBarPMTTimeBottom->Draw();
    c33->cd(6);
    hAnaBarXZPMT->Draw("COLZ");

    plotDetector(hAnaBarXZPMT);
    plotSinglePoints(hAnaBarXZPMT);

    c33->DrawClone();
    c33->Print("plots/c33.pdf");

    return c33;

}

TCanvas* plotC34(){

    auto hAnaBarXZPMT = v[0].Histo2D({"h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200},"anaBarXPMT","anaBarZPMT");

    TCanvas* c34 = new TCanvas("c34","c34",800,800);
    c34->Divide(1,1,0.01,0.01,0);

    c34->cd(1);
    hAnaBarXZPMT->Draw("COLZ");

    plotDetector(hAnaBarXZPMT);
    plotSinglePoints(hAnaBarXZPMT);

    c34->DrawClone();
    c34->Print("plots/c34.pdf");

    return c34;

}

TCanvas* plotC3(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hDetectorNhits = v[0].Histo1D("Detector_Nhits");
	auto hDetectorPdg = v[0].Histo1D("anaBarPDG");
	auto hDetectorID = v[0].Histo1D("anaBarID");
	auto hFingerPdg = v[0].Histo1D("fingerPDG");
	auto hFingerID = v[0].Histo1D("fingerID");
	auto hPMTID = v[0].Histo1D("PMT_id");
	auto hAnaBarPMTID = v[0].Histo1D("anaBarPMTID");
	auto hFingerPMTID = v[0].Histo1D("fingerPMTID");
	auto hAnaBarPMTTime = v[0].Histo1D("anaBarPMTTime");

	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	c3->Divide(3,3,0.01,0.01,0);

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
	c3->cd(7);
	hFingerPMTID->Draw();
	c3->cd(8);
	hAnaBarPMTID->Draw();
	c3->cd(9);
	hAnaBarPMTTime->Draw();

	c3->DrawClone();
	c3->Print("plots/c3RA.pdf");

	return c3;

}


TCanvas* plotC4(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFingerEd = v[0].Histo1D("fingerEd");
	auto hFingerPMTNphot = v[0].Histo1D("fingerPMTNPhotons");
	auto hAnaBarPMTNphot = v[0].Histo1D("anaBarPMTNPhotons");
	auto hAnaBarEd = v[0].Histo1D("anaBarEd");

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

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hAnaBarX = v[0].Histo1D("anaBarXVec");
	auto hAnaBarY = v[0].Histo1D("anaBarYVec");
	auto hAnaBarZ = v[0].Histo1D("anaBarZVec");
	auto hAnaBarT = v[0].Histo1D("anaBarTVec");

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

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hE1vsE2 = v[0].Histo2D({"h2", "E1 vs E2", 100, 0.01, 10.0, 100, 0.01, 30.0},"fingerEd","anaBarEdTotal");

	TCanvas *c6 = new TCanvas("c6","c6",800,800);
	c6->Divide(1,1,0.01,0.01,0);

	c6->cd(1);
	hE1vsE2->Draw("COLZ");

	c6->DrawClone();
	c6->Print("plots/c6RA.pdf");

	return c6;

}

TCanvas* plotC7(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFinger_Edep_vs_Nphot = v[0].Filter("trigger2").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot = v[0].Filter("trigger2").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1 = v[0].Filter("trigger2").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

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

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFinger_Edep_vs_NphotCut = v[0].Filter("trigger3").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_NphotCut = v[0].Filter("trigger3").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1Cut = v[0].Filter("trigger3").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

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

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hAnaBarMult = v[0].Histo1D("imult");

	TCanvas *c11 = new TCanvas("c11", "c11", 800,800);
	c11->Divide(1,1, 0.01, 0.01, 0);

	c11->cd(1);
	hAnaBarMult->Draw();

	c11->DrawClone();
	c11->Print("plots/c11RA.pdf");

	return c11;

}


TCanvas* plotC12(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimPx = v[0].Histo1D("fPx");
	auto hPrimPy = v[0].Histo1D("fPy");
	auto hPrimPz = v[0].Histo1D("fPz");

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

TCanvas* plotC13() {


	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "G4SBS Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0},"anaBarXVec","fPx");
	auto hPz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "G4SBS Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0},"anaBarZVec","fPz");
	auto hz_vs_x = v[0].Filter("trigger2").Histo2D({"h35", "CDet z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0},"anaBarXVec","anaBarZVec");
	auto hPrimXZ = v[0].Histo2D({"h99", "G4SBS z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");



	TCanvas *c13 = new TCanvas("c13","c13",800,800);
	c13->Divide(2,2,0.01,0.01,0);

	c13->cd(1);
	hPx_vs_x->Draw("COLZ");
	c13->cd(2);
	hPz_vs_z->Draw("COLZ");
	c13->cd(3);
	hz_vs_x->Draw("COLZ");
	c13->cd(4);
	hPrimXZ->Draw("COLZ");
    plotDetector(hPrimXZ);

	c13->DrawClone();
	c13->Print("plots/c13.pdf");

	return c13;

}


TCanvas* plotC14(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimX = v[0].Histo1D("Prim_X");
	auto hPrimY = v[0].Histo1D("Prim_Y");
	auto hPrimZ = v[0].Histo1D("Prim_Z");
	auto hPrimXZ = v[0].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");
	
	TCanvas *c14 = new TCanvas("c14","c14",800,800);
	c14->Divide(2,2,0.01,0.01,0);

	c14->cd(1);
	hPrimZ->Draw();
	c14->cd(2);
	hPrimY->Draw();
	c14->cd(3);
	hPrimZ->Draw();
	c14->cd(4);
	hPrimXZ->Draw("COLZ");
        plotDetector(hPrimXZ);

	c14->DrawClone();
	c14->Print("plots/c14.pdf");

	return c14;

}

TCanvas* plotC15() {


	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "Px vs x", 100, -80.0, 80.0, 100, -1.0, 1.0},"anaBarXPMT","anaBarPXPMT");
	auto hPz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "Pz vs z", 100, -240.0, 240.0, 100, -1.0, 1.0},"anaBarZPMT","anaBarPZPMT");
	auto hz_vs_x = v[0].Filter("trigger2").Histo2D({"h35", "z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0},"anaBarXPMT","anaBarZPMT");
	auto hPrimXZ = v[0].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");



	TCanvas *c15 = new TCanvas("c15","c15",800,800);
	c15->Divide(2,2,0.01,0.01,0);

	c15->cd(1);
	hPx_vs_x->Draw("COLZ");
	c15->cd(2);
	hPz_vs_z->Draw("COLZ");
	c15->cd(3);
	hz_vs_x->Draw("COLZ");
	c15->cd(4);
	hPrimXZ->Draw("COLZ");
        plotDetector(hPrimXZ);

	c15->DrawClone();
	c15->Print("plots/c15.pdf");

	return c15;

}

TCanvas* plotC16() {


	auto hx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "X vs XPMT", 100, -80.0, 80.0, 100, -80.0, 80.0},"anaBarXPMT","Prim_X");
	auto hz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "Z vs ZPMT", 100, -240.0, 240.0, 100, -240.0, 240.0},"anaBarZPMT","Prim_Z");
	auto hPx_vs_Px = v[0].Filter("trigger2").Histo2D({"h35", "Px vs PxPMT", 100, -1.0, 1.0, 100, -800.0, 800.0},"anaBarPXPMT","fPx");
	auto hPz_vs_Pz = v[0].Filter("trigger2").Histo2D({"h99", "Pz vs PzPMT", 100, -0.6, 0.6, 100, -2400.0, 2400.0},"anaBarPZPMT","fPz");

	TCanvas *c16 = new TCanvas("c16","c16",800,800);
	c16->Divide(2,2,0.01,0.01,0);

	c16->cd(1);
	hx_vs_x->Draw("COLZ");
	c16->cd(2);
	hz_vs_z->Draw("COLZ");
	c16->cd(3);
	hPx_vs_Px->Draw("COLZ");
	c16->cd(4);
	hPz_vs_Pz->Draw("COLZ");

	c16->DrawClone();
	c16->Print("plots/c16.pdf");

	return c16;

}
