// SBS a Geant-4 Based Model of Hall-A Experiments with 11 GeV
// J.R.M Annand, University of Glasgow
// Class PolNucleonRotate
// Azimuthal modulation of polarised nucleon scattering
// 20/06/09 Adapt Derek Glazier's class DGPolHadronElasticPhysics
// 09/04/10 JRMA check updates to equivalent G4 class.
// 17/04/10 JRMA New class separated from hadronic processes
// 15/11/10 JRMA Bug fix theta in A_y, add Dubna-neutron
// 20/11/10 JRMA Save polarisation
//  4/02/12 JRMA Extend range of models and different p and 12C models
// 29/01/16 JRMA Fix bug polar scattering angle
// 03/02/16 JRMA Extend checks of type of scattering
 
#ifndef PolNucleonRotate_h
#define PolNucleonRotate_h 1

#include "G4HadronicProcess.hh"
#include "G4SystemOfUnits.hh"
enum { EppEl, EnpEl, EnpCEEl, EppInel, EnpInel, EnpCEInel, EnnInel };
enum { EAyConst, EAyLadygin, EAySpinka, EAyMcNaughton, EAyGlister, EAyAzhgirey,
       EAyCE};
//enum { Epp, Enp, EpC, EnC };
static const G4double mp = 0.938272;  // proton mass
static const G4double mn = 0.939565;  // neutron mass


class PolNucleonRotate
{
public:
    
  PolNucleonRotate( G4int ); 
  //  virtual ~PolNucleonRotate();
private:
  G4double fppElpm;
  G4double fnpElpm;
  G4double fnpCEElpm;
  G4double fppInelpm;
  G4double fpnInelpm;
  G4double fnpInelpm;
  G4double fnnInelpm;
  G4int fAyModel[EnnInel+1];
  G4double fAyFactor[EnnInel+1];
  G4LorentzVector fP0;      // incident 4-momentum
  G4LorentzVector fP1;      // recoil 4-momentum
  G4LorentzVector fP2;      // recoil 4-momentum
  G4double fp0;
  G4double fp1;
  G4double fp2;
  G4double ftmin;
  G4int fPDG0;              // ID of primary particle
  G4int fPDG1;              // ID of secondary particle
  G4int fPDG2;              // ID of secondary particle
  G4int fChannel;
  G4bool fIsEl;             // elastic?


public:
  G4double  GetPolarisedRotation(const G4DynamicParticle*,
  				 const G4DynamicParticle*,
  				 const G4DynamicParticle* = NULL,
  				 G4bool = false );
  G4double Ay();
  // particular parametrisations of Ay
  G4double AyLadygin( G4double, G4double );             // free np,pp
  G4double AySpinka( G4double, G4double );              // free pp
  G4double AyMcNaughton( G4double, G4double );          // p-12C
  G4double AyGlister( G4double, G4double );             // p-12C
  G4double AyAzhgirey( G4double, G4double );            // p-12C, p-CH
  G4double AyCE( G4double );                            // high energy CE np
  void SetAyModel(G4int, G4int, G4double);              // set Ay models
  G4double Get_t( G4double, G4double, G4double, G4double);
  G4int fVerbose;

  G4ThreeVector* PolNN;
  G4double gAy;
  G4int gChannel;

};

inline G4double PolNucleonRotate::Ay()
{
  if(fChannel < 0)
    return 0.0;
  G4int model = fAyModel[fChannel];
  G4double ay;
  G4double t = (fP1 - fP0).mag2()/(GeV*GeV);     // 4 momentum transfer
  G4double th = fP0.angle(fP1.vect());           // scattering angle
  G4double pt = (fP0.vect()).mag() * sin(th);    // notional transverse mom
  switch( model ){
  case EAyConst:
  default:
    ay = 1.0;
    break;
  case EAyLadygin:
    ay = AyLadygin(fp0,t);
    break;
  case EAySpinka:
    ay = AySpinka(fp0,t);
    break;
  case EAyMcNaughton:
    ay = AyMcNaughton(fp0,th);
    break;
  case EAyGlister:
    ay = AyGlister(fp0,th);
    break;
  case EAyAzhgirey:
    ay = AyAzhgirey(fp0,th);
    break;
  case EAyCE:
    ay = AyCE(pt);
    break;
  }
  ay *= fAyFactor[fChannel];
  return ay;
}


#endif
