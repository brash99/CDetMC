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
                //cout << "j = " << j << " and Detector_x[j]: " << Detector_x[j] << endl;
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

std::vector<float> getAnaBarPMTTimeOverThreshold(bool trigger, int* PMT_Nphotons, float* PMT_TimeOverThreshold) {

    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttimeoverthreshold[NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS];
    //std::cout << "--------------------" << std::endl;
    if (trigger) {
        for (Int_t icount = AnaBar_PMT_Offset;icount<AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;icount++){
            if (PMT_Nphotons[icount]>Photon_min_cut) {
                //std::cout << "getAnaBarPMTTime: " << icount << " " << PMT_Time[icount] << " " << PMT_Nphotons[icount] << std::endl;
                pmttimeoverthreshold[icount] = PMT_TimeOverThreshold[icount];
                v.push_back(pmttimeoverthreshold[icount]);
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
