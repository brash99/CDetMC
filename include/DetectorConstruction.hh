#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"

//---------------------------------------------------------------------------

class WLSMaterials;
class G4Material;

class G4VPhysicalVolume;
class DetectorSD;
class PMTSD;
class DetectorMessenger;

//---------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* BarConstruction(G4Material* Pscint,G4LogicalVolume* expHall_log,G4Material* Glass);
  G4VPhysicalVolume* FingerConstruction(G4Material* Pscint,G4LogicalVolume* expHall_log,G4Material* Glass);

  void UpdateGeometry();

  inline G4VPhysicalVolume* GetExpHall()    { return fExpHall;  };
  inline G4VPhysicalVolume* GetDetVol()    { return fDetVol;  };
  inline G4VPhysicalVolume* GetDet15Vol()    { return fDet15Vol;  };
  inline G4VPhysicalVolume* GetFingerVol()    { return FingerCounter;  };
  inline G4VPhysicalVolume* GetAnaBarVol()    { return AnaBar;  };
  inline G4VPhysicalVolume* GetCladVol()    { return physiClad;  };
  inline G4VPhysicalVolume* GetWLSfiberVol()    { return physiWLSfiber;  };
  inline DetectorSD*        GetDetSD()      { return fDetSD;    };
  inline PMTSD*             GetPMTSD()      { return fPMTSD;    };

  inline void SetTumourOn     ( G4int    tumon )   { fTumourOn     = tumon; }
  inline void SetTumourRadius ( G4double radius )  { fTumourRadius = radius; }
  inline void SetTumourHeight ( G4double height )  { fTumourHeight = height; }
  inline void SetAnaBarXpos   ( G4double AnaBarXpos )    { fAnaBarXpos   = AnaBarXpos; }
  inline void SetNumberOfLayers   ( G4int NumberOfLayers )    { fNumberOfLayers   = NumberOfLayers; }
  inline void SetAnaBarLength   ( G4double AnaBarLength )    { fAnaBarLength   = AnaBarLength; }
  inline void SetAnaBarWidth   ( G4double AnaBarWidth )    { fAnaBarWidth   = AnaBarWidth; }
  inline void SetAnaBarThickness   ( G4double AnaBarThickness )    { fAnaBarThickness   = AnaBarThickness; }
  inline void SetFingerLength   ( G4double FingerLength )    { fFingerLength   = FingerLength; }
  inline void SetFingerWidth   ( G4double FingerWidth )    { fFingerWidth   = FingerWidth; }
  inline void SetFingerThickness   ( G4double FingerThickness )    { fFingerThickness   = FingerThickness; }
  inline void SetFingerZoffset   ( G4double FingerZoffset )    { fFingerZoffset   = FingerZoffset; }
  inline void SetFingerYoffset   ( G4double FingerYoffset )    { fFingerYoffset   = FingerYoffset; }
  inline void SetFibreDiameter   ( G4double FibreDiameter )    { fFibreDiameter  = FibreDiameter; }
  inline void SetFibreLength   ( G4double FibreLength )    { fFibreLength  = FibreLength; }
  inline void SetHoleDiameter   ( G4double HoleDiameter )    { fHoleDiameter  = HoleDiameter; }
  inline void SetHoleLength   ( G4double HoleLength )    { fHoleLength  = HoleLength; }
  inline void SetCladdingDiameter   ( G4double CladdingDiameter )    { fCladdingDiameter  = CladdingDiameter; }
  inline void SetCladdingLength   ( G4double CladdingLength )    { fCladdingLength  = CladdingLength; }
  inline void SetPhotoCathodeDiameter   ( G4double PhotoCathodeDiameter )    { fPhotoCathodeDiameter  = PhotoCathodeDiameter; }
  inline void SetPhotoCathodeThickness   ( G4double PhotoCathodeThickness)    { fPhotoCathodeThickness  = PhotoCathodeThickness; }
  inline void SetMirrorThickness   ( G4double MirrorThickness)    { fMirrorThickness  = MirrorThickness; }
  inline void SetMylarThickness   ( G4double MylarThickness)    { fMylarThickness  = MylarThickness; }

  G4Material* FindMaterial(G4String);
  
  private:

  WLSMaterials* fMaterials;

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDetVol;
  G4VPhysicalVolume* fDet15Vol;
  G4VPhysicalVolume* FingerCounter;
  G4VPhysicalVolume* AnaBar;
  G4VPhysicalVolume* MylarTop;
  G4VPhysicalVolume* MylarBottom;
  G4VPhysicalVolume* MylarSideFront;
  G4VPhysicalVolume* MylarSideBack;
  G4VPhysicalVolume* MylarFingerFront;
  G4VPhysicalVolume* MylarFingerBack;
  G4VPhysicalVolume* MylarFingerSide1;
  G4VPhysicalVolume* MylarFingerSide2;
  G4VPhysicalVolume* Mirror;
  G4VPhysicalVolume* physiClad;
  G4VPhysicalVolume* physiWLSfiber;

  DetectorSD*        fDetSD;
  PMTSD*             fPMTSD;

  G4int              fTumourOn;
  G4double           fTumourRadius;
  G4double           fTumourHeight;
  G4double           fAnaBarXpos;

  G4int		     fNumberOfLayers;
  G4double           fAnaBarLength;
  G4double           fAnaBarWidth;
  G4double           fAnaBarThickness;
  G4double           fFingerLength;
  G4double           fFingerWidth;
  G4double           fFingerThickness;
  G4double           fFingerZoffset;
  G4double           fFingerYoffset;
  G4double           fFibreDiameter;
  G4double           fFibreLength;
  G4double           fHoleDiameter;
  G4double           fHoleLength;
  G4double           fCladdingDiameter;
  G4double           fCladdingLength;
  G4double           fPhotoCathodeDiameter;
  G4double           fPhotoCathodeThickness;
  G4double           fMirrorThickness;
  G4double           fMylarThickness;

};
#endif

//---------------------------------------------------------------------------

