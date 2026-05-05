#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4Step;

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void SetGeneratedKinematics(G4double thetaEE_deg,
                                G4double thetaEm_deg,
                                G4double kinEm_MeV,
                                G4double thetaEp_deg,
                                G4double kinEp_MeV,
                                const G4ThreeVector& dirEm,
                                const G4ThreeVector& dirEp,
                                G4int componentID);

    void RecordDetectorStep(const G4Step* step,
                            G4int detID,
                            G4int volumeID,
                            G4bool isBoundary);

private:
    void ResetEvent();

    G4double OpeningAngleDeg(const G4ThreeVector& a,
                             const G4ThreeVector& b) const;

private:
    G4int fEventID;

    // Generator-level truth.
    G4double fThetaEEGenDeg;
    G4double fThetaEmGenDeg;
    G4double fThetaEpGenDeg;
    G4double fKinEmGenMeV;
    G4double fKinEpGenMeV;

    G4ThreeVector fDirEmGen;
    G4ThreeVector fDirEpGen;

    // Event component: 0 = IPC-like smooth background, 1 = X17-like signal.
    G4int fComponentID;

    // Ideal reconstruction from generated directions.
    G4double fThetaEEIdealDeg;
    G4double fThetaEEIdealResDeg;

    // Detector-level information.
    G4int fHasEmSi;
    G4int fHasEpSi;
    G4int fHasEmScint;
    G4int fHasEpScint;

    G4int fDetEmSi;
    G4int fDetEpSi;
    G4int fDetEmScint;
    G4int fDetEpScint;

    G4ThreeVector fPosEmSi;
    G4ThreeVector fPosEpSi;
    G4ThreeVector fPosEmScint;
    G4ThreeVector fPosEpScint;

    G4double fThetaEEDetectorRawDeg;

    // Scintillator energy deposits.
    G4double fEdepEmScintMeV;
    G4double fEdepEpScintMeV;
    G4double fEdepAllScintMeV;

    G4int fNHitRows;
    G4int fNSiHitRows;
    G4int fNScintHitRows;
};

#endif
