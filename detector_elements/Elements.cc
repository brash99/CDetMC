/*
 * Elements.cpp
 *
 *  Created on: Nov 17, 2020
 *      Author: petrstepanov
 */

#include "Elements.hh"
// #include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4NistManager.hh"
#include <stdlib.h>

Elements* Elements::fgInstance = NULL;

Elements::Elements() : fElementsList{} {
  // Lead
//  G4Isotope* Pb204 = new G4Isotope("204Pb", 82, 204, 203.97*g/mole);
//  G4Isotope* Pb206 = new G4Isotope("204Pb", 82, 206, 205.97*g/mole);
//  G4Isotope* Pb207 = new G4Isotope("204Pb", 82, 207, 206.98*g/mole);
//  G4Isotope* Pb208 = new G4Isotope("204Pb", 82, 208, 207.98*g/mole);
//  G4Element* Pb = new G4Element("Lead","Pb", 4); // nIsotopes
//  Pb->AddIsotope(Pb204, 1.4000*perCent);
//  Pb->AddIsotope(Pb206, 24.100*perCent);
//  Pb->AddIsotope(Pb207, 22.100*perCent);
//  Pb->AddIsotope(Pb208, 52.400*perCent);
//  fElementsList.push_back(Pb);
//
//  // Tungsten
//  G4Isotope *W180 = new G4Isotope("180W", 74, 180, 179.95*g/mole);
//  G4Isotope *W182 = new G4Isotope("182W", 74, 182, 181.95*g/mole);
//  G4Isotope *W183 = new G4Isotope("183W", 74, 183, 182.95*g/mole);
//  G4Isotope *W184 = new G4Isotope("184W", 74, 184, 183.95*g/mole);
//  G4Isotope *W186 = new G4Isotope("186W", 74, 186, 185.95*g/mole);
//  G4Element* W  = new G4Element("Tungsten","W" , 5); //Introducimos el Wolframio
//  W->AddIsotope(W180, 0.1200*perCent);
//  W->AddIsotope(W182, 26.500*perCent);
//  W->AddIsotope(W183, 14.310*perCent);
//  W->AddIsotope(W184, 30.640*perCent);
//  W->AddIsotope(W186, 28.430*perCent);
//  fElementsList.push_back(W);

  // Trying to use the NIST database

  // Lead
  G4Element* elPb = G4NistManager::Instance()->FindOrBuildElement("Pb");
  fElementsList.push_back(elPb);

  // Tungsten
  G4Element* elW = G4NistManager::Instance()->FindOrBuildElement("W");
  fElementsList.push_back(elW);

  // Nitrogen
  G4Element* elN = G4NistManager::Instance()->FindOrBuildElement("N");
  fElementsList.push_back(elN);

  // Oxygen
  G4Element* elO = G4NistManager::Instance()->FindOrBuildElement("O");
  fElementsList.push_back(elO);

  // Aluminum
  G4Element* elAl = G4NistManager::Instance()->FindOrBuildElement("Al");
  fElementsList.push_back(elAl);

  // Carbon
  G4Element* elC = G4NistManager::Instance()->FindOrBuildElement("C");
  fElementsList.push_back(elC);

  // Hydrogen
  G4Element* elH = G4NistManager::Instance()->FindOrBuildElement("H");
  fElementsList.push_back(elH);

  // Nickel
  G4Element* elNi = G4NistManager::Instance()->FindOrBuildElement("Ni");
  fElementsList.push_back(elNi);

  // Iron
  G4Element* elFe = G4NistManager::Instance()->FindOrBuildElement("Fe");
  fElementsList.push_back(elFe);

  // Copper
  G4Element* elCu = G4NistManager::Instance()->FindOrBuildElement("Cu");
  fElementsList.push_back(elCu);

  // Chromium
  G4Element* elCr = G4NistManager::Instance()->FindOrBuildElement("Cr");
  fElementsList.push_back(elCr);

  // Barium
  G4Element* elBa = G4NistManager::Instance()->FindOrBuildElement("Ba");
  fElementsList.push_back(elBa);

  // Gadolinium
  G4Element* elGd = G4NistManager::Instance()->FindOrBuildElement("Gd");
  fElementsList.push_back(elGd);

  // Germanium
  G4Element* elGe = G4NistManager::Instance()->FindOrBuildElement("Ge");
  fElementsList.push_back(elGe);

  // Silicon
  G4Element* elSi = G4NistManager::Instance()->FindOrBuildElement("Si");
  fElementsList.push_back(elSi);
}

Elements::~Elements() {
}

G4Element* Elements::getElement(const char* symbol){
  for (G4int i=0; i<(int)fElementsList.size(); i++){
    if (fElementsList[i]->GetSymbol() == symbol){
      return fElementsList[i];
    }
  }
  G4cout << "ERROR: Elements::getElement element \"" << symbol << "\" not found";
  exit(1);
  return NULL;
}

Elements* Elements::getInstance() {
  static Elements instance;
  if (fgInstance == nullptr){
    fgInstance = &instance;
  }
  return fgInstance;
}

// PS: when having a static std::vector member
// std::vector<G4Element*> Elements::fElementsList = {};
