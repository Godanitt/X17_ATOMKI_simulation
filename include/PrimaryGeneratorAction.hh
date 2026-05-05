#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <vector>

class EventAction;
class G4Event;
class G4ParticleGun;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(EventAction* eventAction,
                           const G4String& inputFile);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    struct PairRow
    {
        G4int dataRow = -1;                  // 1-based line number in data/data_pair_creation.txt
        G4double thetaEE_deg = 0.0;          // column 1: opening angle e-/e+
        G4double thetaEe_deg = 0.0;          // column 2: electron e- polar angle
        G4double energyEe_MeV = 0.0;         // column 3: electron e- kinetic energy
        G4double thetaEp_deg = 0.0;          // column 4: positron e+ polar angle
        G4double energyEp_MeV = 0.0;         // column 5: positron e+ kinetic energy
    };

    void LoadInputTable(const G4String& inputFile);
    G4ThreeVector DirectionFromPolarAzimuth(G4double theta, G4double phi) const;
    G4double WrapPhiDeg(G4double phiDeg) const;

private:
    G4ParticleGun* fGun = nullptr;
    EventAction* fEventAction = nullptr;
    std::vector<PairRow> fRows;
};

#endif
