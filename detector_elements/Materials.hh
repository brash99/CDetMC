/*
 * Materials.h
 *
 *  Created on: Nov 16, 2020
 *      Author: petrstepanov
 */

#ifndef SRC_MATERIALS_HH_
#define SRC_MATERIALS_HH_

#include "G4Material.hh"

#include <G4MaterialPropertyVector.hh>
#include <vector>

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wshadow"
//#include <TGraph.h>
//#pragma GCC diagnostic pop

enum class BaseMaterialType {
  PWO, SciGlass, SciGlassLight, CSGlass, Dummy
};

class Materials {
public:
  static Materials* getInstance();
  virtual ~Materials();

  G4Material* getMaterial(const char*);

  void printMaterialProperties(const char*);
  void printMaterialProperties(G4Material *material);

  // void setCrystalLength(G4double crystalLength);

protected:
  Materials();

private:
  static Materials *fgInstance;

  void generateEmissionMacro(G4Material *material);

  // void initMaterials();
  // void initCrystalOpticalProperties();

  std::vector<G4Material*> fMaterialsList;
  // List of defined materials
  // static std::vector<G4Material*> fMaterialsList;

//  std::vector<G4double> fOpticalPhotonWavelength;
//  std::vector<G4double> fOpticalPhotonEnergy;

  G4MaterialPropertyVector* calcAbsorptionLength(G4MaterialPropertiesTable *mpt,
      G4MaterialPropertyVector *refractiveIndex,
      G4MaterialPropertyVector *measuredLongitudinalTransmittance,
      G4double sampleLength, G4int nPoints = 50);

//  void printLine(const char* heading, std::vector<G4double> values, G4double unit = 0.);
//  void printHeader();

  G4MaterialPropertyVector* nmToMPV(std::vector<G4double> digitizedArrayNm);
  G4MaterialPropertyVector* evToMPV(std::vector<G4double> digitizedArrayNm);

  G4double wavelengthNmToEnergy(G4double wavelength);
  G4double energyToWavelangthNm(G4double energy);

  const G4int TAB_COLUMN_1 = 25;
  const G4int TAB_COLUMN = 10;

  G4Material* getBaseMaterial(const char *newName, BaseMaterialType bmt,
      G4double overrideDensity = 0);
  G4MaterialPropertiesTable* getBaseMPT(BaseMaterialType bmt);

protected:
//  TGraph* drawMPV(G4MaterialPropertyVector *mpv);
//  void saveMaterial(G4Material *material);
};

#endif /* SRC_MATERIALS_HH_ */
