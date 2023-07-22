#ifndef PMTHit_h
#define PMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class PMTHit : public G4VHit
{
public:
  
  PMTHit();
  ~PMTHit();
  PMTHit(const PMTHit &right);
  const PMTHit& operator=(const PMTHit &right);

  G4int operator==(const PMTHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

  inline void IncPhotonCount(){photons++;}
  inline G4int GetPhotonCount(){return photons;}

  inline void SetPMTNumber(G4int n) { pmtNumber = n; }
  inline G4int GetPMTNumber() { return pmtNumber; }

  inline void SetPMTTime(G4double t) { pmtTime = t; }
  inline G4double GetPMTTime() { return pmtTime; }

  //inline void SetPMTKineticEnergy(G4int photonhits, G4double energy) { pmtKineticEnergy[photonhits] = energy; }
  //inline G4double GetPMTKineticEnergy(G4int photonhits) { return pmtKineticEnergy[photonhits]; }

private:
  G4int pmtNumber;
  G4int photons;
  G4double pmtTime;
  //G4double pmtKineticEnergy[50000];
};

typedef G4THitsCollection<PMTHit> PMTHitsCollection;

extern G4Allocator<PMTHit> PMTHitAllocator;

inline void* PMTHit::operator new(size_t){
  void *aHit;
  aHit = (void *) PMTHitAllocator.MallocSingle();
  return aHit;
}

inline void PMTHit::operator delete(void *aHit){
  PMTHitAllocator.FreeSingle((PMTHit*) aHit);
}

#endif


