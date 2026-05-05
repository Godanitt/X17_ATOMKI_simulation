#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <vector>

class G4ParticleGun;
class G4Event;
class EventAction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(EventAction* eventAction,
                           const G4String& filename = "data/data_pair_creation.txt");

    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    enum EventComponent
    {
        kIPCBackground = 0,
        kX17Signal = 1,
        kAccidentalReserved = 2
    };

    struct PairEvent
    {
        G4double thetaEE_deg;
        G4double thetaElectron_deg;
        G4double kineticElectron_MeV;
        G4double thetaPositron_deg;
        G4double kineticPositron_MeV;
    };

    void LoadTable(const G4String& filename);
    void ConfigureFromEnvironment();

    PairEvent SampleSignalPair() const;
    PairEvent SampleIPCBackgroundPair() const;

    G4double SampleTruncatedExponentialDeg(G4double thetaMinDeg,
                                           G4double thetaMaxDeg,
                                           G4double slopePerDeg) const;

    G4ThreeVector DirectionFromThetaPhi(G4double theta,
                                        G4double phi) const;

private:
    G4ParticleGun* fGun;
    EventAction* fEventAction;
    std::vector<PairEvent> fTable;

    // Generator mixture. IPC is deliberately dominant by default because it is
    // the smooth physical background under the X17-like large-angle structure.
    G4double fSignalFraction;
    G4double fIPCThetaMinDeg;
    G4double fIPCThetaMaxDeg;
    G4double fIPCThetaSlopePerDeg;
    G4double fIPCEnergyAsymmetryMax;
    G4double fIPCEnergyScaleSigma;
};

#endif
