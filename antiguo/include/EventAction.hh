#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
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

    // Direct copy of the selected line of data/data_pair_creation.txt:
    //   col1 thetaEE, col2 thetaEe, col3 energyEe, col4 thetaEp, col5 energyEp.
    void SetGeneratedKinematics(G4int dataRow,
                                G4double thetaEEInputDeg,
                                G4double thetaEeInputDeg,
                                G4double thetaEpInputDeg,
                                G4double energyEeMeV,
                                G4double energyEpMeV,
                                G4double phi0Deg,
                                G4double phiEeDeg,
                                G4double phiEpDeg,
                                const G4ThreeVector& dirEe,
                                const G4ThreeVector& dirEp);

    void RecordHit(const G4Step* step,
                   G4int detID,
                   G4int volumeID);

private:
    void Reset();
    G4double OpeningAngleDeg(const G4ThreeVector& a,
                             const G4ThreeVector& b) const;

private:
    G4int fEventID = -1;
    G4int fDataRow = -1;

    G4double fThetaEEInputDeg = -999.0;
    G4double fThetaEeInputDeg = -999.0;
    G4double fThetaEpInputDeg = -999.0;
    G4double fEnergyEeMeV = -999.0;
    G4double fEnergyEpMeV = -999.0;
    G4double fPhi0Deg = -999.0;
    G4double fPhiEeDeg = -999.0;
    G4double fPhiEpDeg = -999.0;

    G4ThreeVector fDirEe;
    G4ThreeVector fDirEp;
};

#endif
