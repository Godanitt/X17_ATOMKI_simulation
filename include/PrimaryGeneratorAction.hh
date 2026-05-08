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
                           const G4String& inputFile,
                           const G4String& generationMode);
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

    enum class Mode { SignalFromTxt, SmoothBackground };

    void LoadInputTable(const G4String& inputFile);
    void GenerateSignalFromTxt(G4Event* event, G4int eventID);
    void GenerateSmoothBackground(G4Event* event, G4int eventID);

    G4ThreeVector DirectionFromPolarAzimuth(G4double theta, G4double phi) const;
    G4ThreeVector RandomIsotropicDirection() const;
    G4ThreeVector DirectionAtOpeningAngle(const G4ThreeVector& axis,
                                          G4double openingAngle,
                                          G4double rotationAngle) const;
    G4double SampleTruncatedExponentialDeg(G4double thetaMaxDeg,
                                           G4double thetaScaleDeg) const;
    G4double WrapPhiDeg(G4double phiDeg) const;

private:
    G4ParticleGun* fGun = nullptr;
    EventAction* fEventAction = nullptr;
    Mode fMode = Mode::SignalFromTxt;
    std::vector<PairRow> fRows;

    // Simple smooth IPC-like background model. It is deliberately not tuned to
    // data: it only provides a broad falling opening-angle template that can be
    // passed through the same geometry and detector-effects chain as the signal.
    G4double fBackgroundTotalKineticMeV = 17.5;
    G4double fBackgroundMinKineticMeV = 0.20;
    G4double fBackgroundThetaScaleDeg = 40.0;
};

#endif
