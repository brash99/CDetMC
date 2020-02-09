// Initial version from David Hamilton
//
// Initial modifications from Kerver and Brash for geometry of finger
// counter + AnaBar
//
// Additional modifications by Hamilton and Brash in order to add 
// optical photon capabilities, need to include:
// G4OpticalSurface, G4LogicalBorderSurface, and G4LogicalSkinSurface
// In order to record PMT hits, also include PMTSD.hh
//
// Additional mods by Brash to add WLS fibre + cladding for CDet simulation

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "DetectorSD.hh"
#include "PMTSD.hh"
#include "WLSMaterials.hh"

#include "G4Material.hh"
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction()
  : fMaterials(NULL)
{
  fDetMessenger = new DetectorMessenger(this);

  fNistManager  = G4NistManager::Instance();

  fTumourRadius = 0.5;
  fTumourHeight = 0.0;
  fAnaBarXpos	= 0.0;

//TODO increase this to 3x14 and 14x14 to model a half module
  fNumberOfLayers = 14;
  
  fMirrorThickness = 0.20;
  fMylarThickness = 0.02;

  fAnaBarLength = 50.0;
  fAnaBarWidth = 4.0;
  fAnaBarThickness = 0.50;
  //fAnaBarThickness = 5.0;

  //fFingerLength = 2.6;
  //fFingerWidth = 4.0;
  //fFingerThickness = 1.7;
  //fFingerZoffset = -1.0;
 //TODO 2 triggers up top and 2 triggers below. 2xwidth? 
  fFingerLength = 10.0;
  fFingerWidth =  fNumberOfLayers*(fAnaBarThickness+2.0*fMylarThickness)+20.0;
  fFingerThickness = 1.0;
  fFingerZoffset = -(fFingerWidth-20.0)/2.0;
  fFingerYoffset = fAnaBarWidth/2.0+fFingerThickness/2.0+1.0;

  fHoleDiameter = 0.19;
  //fHoleDiameter = 1.6;
  fHoleLength = fAnaBarLength;

  fCladdingDiameter = 0.18;
  //fCladdingDiameter = 1.8;
  fCladdingLength = 60.0;
  
  fFibreDiameter = 0.16;
  //fFibreDiameter = 1.6;
  fFibreLength = fCladdingLength;

  fPhotoCathodeDiameter = 2.54;
  fPhotoCathodeThickness = 0.30;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command = "/control/execute macros/DetectorSetup.mac";
  UI->ApplyCommand(command);

}

//---------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() 
{
  if (fMaterials)	delete fMaterials;
}

//---------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //---------------------------------------------------------------------------
  // Define Materials and Optical Properties
  //---------------------------------------------------------------------------

  fMaterials = WLSMaterials::GetInstance();

  G4Element* H  = new G4Element("H",  "H",  1.,  1.01*g/mole);
  G4Element* B  = new G4Element("B",  "B",  5.,  10.8110*g/mole);
  G4Element* C  = new G4Element("C",  "C",  6.,  12.01*g/mole);
  G4Element* N  = new G4Element("N",  "N",  7.,  14.01*g/mole);
  G4Element* O  = new G4Element("O",  "O",  8.,  16.00*g/mole);
  G4Element* Na = new G4Element("Na", "Na", 11., 22.9898*g/mole);
  G4Element* Si = new G4Element("Si", "Si", 14., 28.0855*g/mole);

  //-----------------------------------------------------
  
  //const G4int Num         = 6;
  //G4double    Energy[Num] = { 2.34*eV , 2.46*eV,  2.58*eV, 2.72*eV,  2.88 *eV, 3.06 *eV};
  const G4int Num         = 12;
  G4double    Energy[Num] = { 3.44*eV, 3.26*eV, 3.1*eV, 3.02*eV, 2.95*eV,
			      2.92*eV, 2.82*eV, 2.76*eV, 2.7*eV, 2.58*eV,
			      2.38*eV, 2.08*eV};

  //-----------------------------------------------------
  // Air
  G4Material* Air = new G4Material("Air", 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  G4double    Air_RIND[Num]         = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. };
  G4double    Air_AbsLength[Num]    = { 4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm,
  					4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm, 4200.*cm };
  G4MaterialPropertiesTable *Air_mt = new G4MaterialPropertiesTable();
  Air_mt->AddProperty("RINDEX",    Energy, Air_RIND,      Num );
  Air_mt->AddProperty("ABSLENGTH", Energy, Air_AbsLength, Num );
  Air->SetMaterialPropertiesTable(Air_mt);

  //-----------------------------------------------------
  // Plastic Scintillator
  //
  // Update properties to model BC-408 ...
  // From https://www-zeuthen.desy.de/lcdet/Feb_05_WS/talks/rd_pcdet_sim.pdf
  //
  G4Material* Pscint = new G4Material("Pscint", 1.03*g/cm3, 2);
  Pscint->AddElement(C, 8);
  Pscint->AddElement(H, 8);

  G4double    Pscint_RIND[Num]      = { 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
  					1.58, 1.58, 1.58, 1.58, 1.58, 1.58 };
  G4double    Pscint_AbsLength[Num] = { 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm,
  					210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm };
  G4double    Pscint_SCINT[Num]     = { 0.04, 0.07, 0.20, 0.49, 0.84, 1.00,
  					0.83, 0.55, 0.40, 0.17, 0.03, 0.005 };

  G4MaterialPropertiesTable *Pscint_mt = new G4MaterialPropertiesTable();
  Pscint_mt->AddProperty("RINDEX",        Energy, Pscint_RIND,      Num );
  Pscint_mt->AddProperty("ABSLENGTH",     Energy, Pscint_AbsLength, Num );
  Pscint_mt->AddProperty("FASTCOMPONENT", Energy, Pscint_SCINT,     Num );
  Pscint_mt->AddProperty("SLOWCOMPONENT", Energy, Pscint_SCINT,     Num );

  Pscint_mt->AddConstProperty("SCINTILLATIONYIELD", 500./MeV ); 
  Pscint_mt->AddConstProperty("RESOLUTIONSCALE" ,   1.0        ); 
  Pscint_mt->AddConstProperty("FASTTIMECONSTANT",   1.0 *ns    );  
  Pscint_mt->AddConstProperty("SLOWTIMECONSTANT",   1.0 *ns    );  
  Pscint_mt->AddConstProperty("YIELDRATIO",         1.0        );

  Pscint->SetMaterialPropertiesTable(Pscint_mt);

  //-----------------------------------------------------
  // PMT Glass

  G4Material* Glass = new G4Material("Glass", 2.55*g/cm3, 4);
  Glass->AddElement(Na, 0.0300);
  Glass->AddElement(B,  0.0400);
  Glass->AddElement(O, 0.5395);
  Glass->AddElement(Si, 0.3905);

  //G4double    Glass_RIND[Num]      = { 1.54367, 1.54332, 1.53440, 1.53114, 1.52770, 1.52562 };
  //G4double    Glass_AbsLength[Num] = { 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm };
  //G4double    Glass_REFL[Num]      = { 0.0909, 0.0898, 0.0864, 0.0867, 0.0853, 0.0846 };
  G4double    Glass_RIND[Num]      = { 1.525, 1.525, 1.525, 1.525, 1.525, 1.528,
  					1.528, 1.528, 1.531, 1.531, 1.534, 1.534 };
  G4double    Glass_AbsLength[Num] = { 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm,
  					420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm, 420.*cm };
  G4double    Glass_REFL[Num]      = { 0.0846, 0.0846, 0.0846, 0.0846, 0.0846, 0.0853,
  					 0.0853, 0.0853, 0.0864, 0.0867, 0.0909, 0.0909 };

  G4MaterialPropertiesTable *Glass_mt = new G4MaterialPropertiesTable();
  Glass_mt->AddProperty("ABSLENGTH",    Energy, Glass_AbsLength, Num);
  Glass_mt->AddProperty("RINDEX",       Energy, Glass_RIND,      Num);
  Glass_mt->AddProperty("REFLECTIVITY", Energy, Glass_REFL,      Num);
  
  Glass->SetMaterialPropertiesTable(Glass_mt);

  //---------------------------------------------------------------------------
  // Create Experimental Hall
  //---------------------------------------------------------------------------

  G4Box* expHall_box           = new G4Box("expHall_box",
					   1.5 *m, 1.5 *m, 1.5 *m );
  
  //G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
  //						     FindMaterial("G4_AIR"),
  //						     "expHall_log", 0, 0, 0);
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,
						     Air,
						     "expHall_log", 0, 0, 0);
  
  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   expHall_log, "expHall", 0, false, 0);

  //----------------- Create Mylar Layers and Sides ----------------------

  G4VSolid* solidMylar = new G4Box("Mylar",fAnaBarLength/2.0*cm,fAnaBarWidth/2.0*cm,fMylarThickness/2.0*cm);
  
  G4LogicalVolume* logicMylar = new G4LogicalVolume(solidMylar,FindMaterial("G4_Al"), "Mylar");
  
  G4VSolid* solidMylarSide = new G4Box("MylarSide",fAnaBarLength/2.0*cm,fMylarThickness/2.0*cm,fNumberOfLayers*(fAnaBarThickness/2.0+fMylarThickness)*cm);
  
  G4LogicalVolume* logicMylarSide = new G4LogicalVolume(solidMylarSide,FindMaterial("G4_Al"), "MylarSide");

  G4OpticalSurface* mylarSurface = new G4OpticalSurface("MylarSurface",glisur, ground, dielectric_metal, 1.0);

  G4MaterialPropertiesTable* mylarSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_mylar[] = {2.00*eV, 3.47*eV};
  const G4int nybins = sizeof(p_mylar)/sizeof(G4double);
  G4double refl_mylar[] = {0.9, 0.9};
  assert (sizeof(refl_mylar) == sizeof(p_mylar));
  G4double effi_mylar[] = {0, 0};
  assert (sizeof(effi_mylar) == sizeof(p_mylar));

  mylarSurfaceProperty->AddProperty("REFLECTIVITY",p_mylar,refl_mylar,nybins);
  mylarSurfaceProperty->AddProperty("EFFICIENCY",p_mylar,effi_mylar,nybins);

  mylarSurface -> SetMaterialPropertiesTable(mylarSurfaceProperty);

  new G4LogicalSkinSurface("MylarSurface",logicMylar,mylarSurface);
  new G4LogicalSkinSurface("MylarSurfaceSide",logicMylarSide,mylarSurface);

  //----------------- Create Mylar surrounding Finger PMT  ----------------------

  G4VSolid* solidMylarFinger = new G4Box("MylarFinger",fFingerLength/2.0*cm,fMylarThickness/2.0*cm,fFingerWidth/2.0*cm);
  
  G4LogicalVolume* logicMylarFinger = new G4LogicalVolume(solidMylarFinger,FindMaterial("G4_Al"), "MylarFinger");
  
  G4VSolid* solidMylarFingerSide = new G4Box("MylarFingerSide",fMylarThickness/2.0*cm,fFingerThickness/2.0*cm,fFingerWidth/2.0*cm);
  
  G4LogicalVolume* logicMylarFingerSide = new G4LogicalVolume(solidMylarFingerSide,FindMaterial("G4_Al"), "MylarFingerSide");

  new G4LogicalSkinSurface("MylarSurfaceFinger",logicMylarFinger,mylarSurface);
  new G4LogicalSkinSurface("MylarSurfaceFingerSide",logicMylarFingerSide,mylarSurface);
  
  G4ThreeVector MylarFinger_pos1(0.0*cm , (fFingerYoffset-fFingerThickness/2.0-fMylarThickness/2.0)*cm , fFingerZoffset*cm);
   	MylarFingerFront    =  new G4PVPlacement(0, MylarFinger_pos1 , logicMylarFinger , "MylarFinger" , expHall_log , false , 82);


  G4ThreeVector MylarFinger_pos2(0.0*cm , (fFingerYoffset+fFingerThickness/2.0+fMylarThickness/2.0)*cm , fFingerZoffset*cm);
   	MylarFingerBack    =  new G4PVPlacement(0, MylarFinger_pos2 , logicMylarFinger , "MylarFinger" , expHall_log , false , 83);

  G4ThreeVector MylarFinger_pos3((fFingerLength/2.0+fMylarThickness/2.0)*cm , fFingerYoffset*cm , fFingerZoffset*cm);
   	MylarFingerSide1    =  new G4PVPlacement(0, MylarFinger_pos3 , logicMylarFingerSide , "MylarFingerSide" , expHall_log , false , 84);

  G4ThreeVector MylarFinger_pos4(-1.0*(fFingerLength/2.0+fMylarThickness/2.0)*cm , fFingerYoffset*cm , fFingerZoffset*cm);
   	MylarFingerSide2    =  new G4PVPlacement(0, MylarFinger_pos4 , logicMylarFingerSide , "MylarFingerSide" , expHall_log , false , 85);


  //---------------------------------------------------------------------------
  // Create Detectors
  //---------------------------------------------------------------------------

   G4Box* fingercounter_solid  = new G4Box("fingercounter_solid", fFingerLength/2.0*cm , fFingerThickness/2.0*cm , fFingerWidth/2.0*cm);
   //G4LogicalVolume* fingercounter_log = new G4LogicalVolume(fingercounter_solid, FindMaterial("Polystyrene"), "fingercounter_log");
   G4LogicalVolume* fingercounter_log = new G4LogicalVolume(fingercounter_solid, Pscint, "fingercounter_log");
   G4ThreeVector fingercounter_pos(0.0*cm , fFingerYoffset*cm , fFingerZoffset*cm);
   FingerCounter=  new G4PVPlacement(0, fingercounter_pos , fingercounter_log , "FingerCounter" , expHall_log , false , 0);

   
   G4VSolid* AnaBar_outer  = new G4Box("AnaBar_solid_outer", fAnaBarLength/2.0*cm , fAnaBarWidth/2.0*cm , fAnaBarThickness/2.0*cm);
   G4VSolid* AnaBar_inner = new G4Tubs("AnaBar_solid_inner", 0.*cm,fHoleDiameter/2.0*cm, fAnaBarLength/2.0*cm, 0.*deg, 360.0*deg);
   
   G4ThreeVector fibre_pos(0.0*cm , 0.0*cm , 0.0*cm);
   G4RotationMatrix* anabar_rm  = new G4RotationMatrix();
   anabar_rm->rotateY(90. *deg);
  
   G4VSolid* AnaBar_solid  = new G4SubtractionSolid("AnaBar_solid", AnaBar_outer, AnaBar_inner, anabar_rm,fibre_pos);

   G4LogicalVolume* AnaBar_log = new G4LogicalVolume(AnaBar_solid, Pscint, "AnaBar_log");
   
   for (G4int iii=0; iii<fNumberOfLayers; iii++){ 
   	G4ThreeVector AnaBar_pos(fAnaBarXpos*cm , 0.0*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*iii*cm);
   	AnaBar      =  new G4PVPlacement(0, AnaBar_pos , AnaBar_log , "AnaBar" , expHall_log , false , iii+1);

   	G4ThreeVector Mylar_pos1(fAnaBarXpos*cm , 0.0*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*iii*cm-fAnaBarThickness/2.0*cm-fMylarThickness/2.0*cm);
   	MylarTop    =  new G4PVPlacement(0, Mylar_pos1 , logicMylar , "Mylar" , expHall_log , false , 50+iii+1);

   	G4ThreeVector Mylar_pos2(fAnaBarXpos*cm , 0.0*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*iii*cm+fAnaBarThickness/2.0*cm+fMylarThickness/2.0*cm);
   	MylarBottom =  new G4PVPlacement(0, Mylar_pos2 , logicMylar , "Mylar" , expHall_log , false , 70+iii+1);
   }

   G4ThreeVector Mylar_pos3(fAnaBarXpos*cm , 1.0*(fAnaBarWidth/2.0+fMylarThickness/2.0)*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*(fNumberOfLayers/2.0-0.5)*cm );
   G4ThreeVector Mylar_pos4(fAnaBarXpos*cm , -1.0*(fAnaBarWidth/2.0+fMylarThickness/2.0)*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*(fNumberOfLayers/2.0-0.5)*cm ); 
   	MylarSideFront =  new G4PVPlacement(0, Mylar_pos3 , logicMylarSide , "MylarSide" , expHall_log , false , 80);
   	MylarSideBack =  new G4PVPlacement(0, Mylar_pos4 , logicMylarSide , "MylarSide" , expHall_log , false , 81);

  //----------------- Create Mirror on non-phototube end ----------------------

  G4VSolid* solidMirror = new G4Box("Mirror",fMirrorThickness/2.0*cm,fAnaBarWidth/2.0*cm,fNumberOfLayers*(fAnaBarThickness/2.0+fMylarThickness)*cm);

  G4LogicalVolume* Mirror_log = new G4LogicalVolume(solidMirror,FindMaterial("G4_Al"), "Mirror");

  G4OpticalSurface* mirrorSurface = new G4OpticalSurface("MirrorSurface",glisur, ground, dielectric_metal, 1.0);

  G4MaterialPropertiesTable* mirrorSurfaceProperty = new G4MaterialPropertiesTable();

  G4double p_mirror[] = {2.00*eV, 3.47*eV};
  const G4int nbins = sizeof(p_mirror)/sizeof(G4double);
  G4double refl_mirror[] = {1.0, 1.0};
  assert (sizeof(refl_mirror) == sizeof(p_mirror));
  G4double effi_mirror[] = {0, 0};
  assert (sizeof(effi_mirror) == sizeof(p_mirror));

  mirrorSurfaceProperty->AddProperty("REFLECTIVITY",p_mirror,refl_mirror,nbins);
  mirrorSurfaceProperty->AddProperty("EFFICIENCY",p_mirror,effi_mirror,nbins);

  mirrorSurface -> SetMaterialPropertiesTable(mirrorSurfaceProperty);

  Mirror = new G4PVPlacement(0,G4ThreeVector(fAnaBarXpos*cm-fAnaBarLength/2.0*cm-fMirrorThickness/2.0*cm, 0.0*cm, -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*(fNumberOfLayers/2.0-0.5)*cm),Mirror_log,"Mirror",expHall_log,false,40);

  new G4LogicalSkinSurface("MirrorSurface",Mirror_log,mirrorSurface);

  //---------------------------------------------------------------------------
  // Create Fibre ... first cladding, and then WLS fibre itself.
  //---------------------------------------------------------------------------
   
  
  G4VSolid* solidClad1;
  
  solidClad1 = new G4Tubs("Clad1",fFibreDiameter/2.0*cm,fCladdingDiameter/2.0*cm,fFibreLength/2.0*cm,0.0*deg,360.0*deg);
  
  G4LogicalVolume* logicClad1 = new G4LogicalVolume(solidClad1, FindMaterial("Pethylene"),"Clad1");

  G4VSolid* solidWLSfiber;

  solidWLSfiber = new G4Tubs("WLSFiber",0.*cm,fFibreDiameter/2.0*cm,fFibreLength/2.0*cm,0.0*deg,360.0*deg);

  G4LogicalVolume* logicWLSfiber = new G4LogicalVolume(solidWLSfiber,
                                                         FindMaterial("PMMA"),
                                                         "WLSFiber");

  logicWLSfiber->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

  for (G4int iii=0; iii<fNumberOfLayers; iii++){ 
  	G4ThreeVector Global_fibre_pos(fAnaBarXpos*cm + (fFibreLength-fAnaBarLength)/2.0*cm , 0.0*cm , -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*iii*cm);

  	physiClad = new G4PVPlacement(anabar_rm,Global_fibre_pos,logicClad1,"Clad1",expHall_log,false,15+2.0*iii);
  	physiWLSfiber = new G4PVPlacement(anabar_rm,
                                       Global_fibre_pos,
                                       logicWLSfiber,
                                       "WLSFiber",
                                       expHall_log,
                                       false,
                                       16+2.0*iii);
  }

  //---------------------------------------------------------------------------
  // Create AnaBar PMT
  //---------------------------------------------------------------------------

  
  G4Tubs* det1_tubs            = new G4Tubs("det1_tubs",
					    0. *mm, fCladdingDiameter/2.0*cm, fPhotoCathodeThickness/2.0*cm,
					    0. *deg, 360. *deg );
  
  G4LogicalVolume* det1_log = new G4LogicalVolume(det1_tubs,
						  Glass,
						  "det1_log", 0, 0, 0);
  for (G4int iii=0; iii<fNumberOfLayers; iii++){
  
  	fDetVol                  = new G4PVPlacement(anabar_rm, G4ThreeVector(fFibreLength/2.0*cm+fAnaBarXpos*cm+fPhotoCathodeThickness/2.0*cm+(fFibreLength/2.0-fAnaBarLength/2.0)*cm, 0., -1.0*(fAnaBarThickness/2.0)*cm-(fAnaBarThickness+2.0*fMylarThickness)*iii*cm),
						det1_log, "det1", expHall_log, false, iii);
  }

  //---------------------------------------------------------------------------
  // Create Finger PMT
  //---------------------------------------------------------------------------
 
       // TODO  Finger(x,y,z) in order to add an extra above the bar and two below
  G4RotationMatrix* finger_rm  = new G4RotationMatrix();
  finger_rm->rotateX(0. *deg);
  
  G4Box* det2_tubs            = new G4Box("det2_tubs",fFingerLength/2.0*cm, fFingerThickness/2.0*cm, fPhotoCathodeThickness/2.0*cm);
  
  G4LogicalVolume* det2_log = new G4LogicalVolume(det2_tubs,
						  Glass,
						  "det2_log", 0, 0, 0);
  
  fDet15Vol                  = new G4PVPlacement(finger_rm, G4ThreeVector(0.0,fFingerYoffset*cm,fFingerZoffset*cm+fFingerWidth/2.0*cm+fPhotoCathodeThickness/20*cm),
						det2_log, "det2", expHall_log, false, 14);

  //---------------------------------------------------------------------------
  // Create Optical Surface
  //---------------------------------------------------------------------------
 
  //G4double PhotoCath_EFF[Num]  = { 0.24, 0.24, 0.24, 0.24, 0.24, 0.24 };
  G4double PhotoCath_EFF[Num]  = { 0.13, 0.13, 0.13, 0.13, 0.13, 0.13 };
  G4double PhotoCath_REFL[Num] = { 0., 0., 0., 0., 0., 0. };
  
  G4MaterialPropertiesTable* PhotoCath_mt = new G4MaterialPropertiesTable();
  PhotoCath_mt->AddProperty("EFFICIENCY",   Energy, PhotoCath_EFF,  Num );
  PhotoCath_mt->AddProperty("REFLECTIVITY", Energy, PhotoCath_REFL, Num );
  
  G4OpticalSurface* PhotoCath_opsurf = new G4OpticalSurface("PhotoCath_opsurf", 
							    glisur, polished, dielectric_metal );
  PhotoCath_opsurf->SetMaterialPropertiesTable(PhotoCath_mt);
  
  new G4LogicalSkinSurface("PhotoCathAnaBar_surf", det1_log, PhotoCath_opsurf );
  new G4LogicalSkinSurface("PhotoCathFinger_surf", det2_log, PhotoCath_opsurf );

  //---------------------------------------------------------------------------
  // Set Step Limits, Sensitive Detector and Visualisation
  //---------------------------------------------------------------------------

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  fDetSD = new DetectorSD("DetSD", 20);
  SDman->AddNewDetector( fDetSD );
  fingercounter_log->SetSensitiveDetector( fDetSD );
  AnaBar_log->SetSensitiveDetector( fDetSD );

  fPMTSD = new PMTSD("PMTSD", 100);
  SDman->AddNewDetector( fPMTSD );
  det1_log->SetSensitiveDetector( fPMTSD );
  det2_log->SetSensitiveDetector( fPMTSD );
  
  G4VisAttributes* blue    = new G4VisAttributes( G4Colour(0.0,0.0,1.0)   );
  G4VisAttributes* yellow  = new G4VisAttributes( G4Colour(0.0,1.0,1.0)   );
  G4VisAttributes* green   = new G4VisAttributes( G4Colour(0.0,1.0,0.0)   );
  G4VisAttributes* red     = new G4VisAttributes( G4Colour(1.0,0.0,0.0)   );
  expHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  fingercounter_log->SetVisAttributes(blue);
  AnaBar_log->SetVisAttributes(blue);
  Mirror_log->SetVisAttributes(yellow);
  logicWLSfiber->SetVisAttributes(green);
  logicMylar->SetVisAttributes(red);
  logicMylarFinger->SetVisAttributes(red);
  logicClad1->SetVisAttributes(yellow);
  logicMylarSide->SetVisAttributes(red);
  logicMylarFingerSide->SetVisAttributes(red);
  det1_log->SetVisAttributes(yellow);
  det2_log->SetVisAttributes(yellow);

  return fExpHall;
}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------

G4Material* DetectorConstruction::FindMaterial(G4String name) {
    std::cout << "Here I am ... rock me like a hurricane! " << name << std::endl;
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}


