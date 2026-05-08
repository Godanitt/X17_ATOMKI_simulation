#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>

EventAction::EventAction()
{
    Reset();
}

void EventAction::Reset()
{
    fEventID = -1;
    fDataRow = -1;

    fThetaEEInputDeg = -999.0;
    fThetaEeInputDeg = -999.0;
    fThetaEpInputDeg = -999.0;
    fEnergyEeMeV = -999.0;
    fEnergyEpMeV = -999.0;
    fPhi0Deg = -999.0;
    fPhiEeDeg = -999.0;
    fPhiEpDeg = -999.0;

    fDirEe = G4ThreeVector(0.0, 0.0, 0.0);
    fDirEp = G4ThreeVector(0.0, 0.0, 0.0);
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    // Do not call Reset() here. In the Geant4 event flow, GeneratePrimaries()
    // may have already copied the txt row into this EventAction object.
    // Resetting here is what produced generated rows full of -999 in older versions.
    fEventID = event ? event->GetEventID() : -1;
}

void EventAction::SetGeneratedKinematics(G4int dataRow,
                                          G4double thetaEEInputDeg,
                                          G4double thetaEeInputDeg,
                                          G4double thetaEpInputDeg,
                                          G4double energyEeMeV,
                                          G4double energyEpMeV,
                                          G4double phi0Deg,
                                          G4double phiEeDeg,
                                          G4double phiEpDeg,
                                          const G4ThreeVector& dirEe,
                                          const G4ThreeVector& dirEp)
{
    fDataRow = dataRow;
    fThetaEEInputDeg = thetaEEInputDeg;
    fThetaEeInputDeg = thetaEeInputDeg;
    fThetaEpInputDeg = thetaEpInputDeg;
    fEnergyEeMeV = energyEeMeV;
    fEnergyEpMeV = energyEpMeV;
    fPhi0Deg = phi0Deg;
    fPhiEeDeg = phiEeDeg;
    fPhiEpDeg = phiEpDeg;
    fDirEe = dirEe.unit();
    fDirEp = dirEp.unit();
}

G4double EventAction::OpeningAngleDeg(const G4ThreeVector& a,
                                      const G4ThreeVector& b) const
{
    if (a.mag2() <= 0.0 || b.mag2() <= 0.0) return -999.0;

    auto cosang = a.unit().dot(b.unit());
    cosang = std::max(-1.0, std::min(1.0, cosang));
    return std::acos(cosang) / deg;
}

void EventAction::RecordHit(const G4Step* step,
                            G4int detID,
                            G4int volumeID)
{
    if (!step) return;

    const auto* track = step->GetTrack();
    const auto* pre = step->GetPreStepPoint();
    if (!track || !pre) return;

    const auto* particle = track->GetParticleDefinition();
    if (!particle) return;

    const auto pdg = particle->GetPDGEncoding();
    if (pdg != 11 && pdg != -11) return;

    const auto pos = pre->GetPosition();
    const auto mom = pre->GetMomentumDirection();

    auto phiHitDeg = pos.phi() / deg;
    if (phiHitDeg < 0.0) phiHitDeg += 360.0;

    const G4int charge = (pdg == 11) ? -1 : +1;
    const G4double thetaGenDeg = (pdg == 11) ? fThetaEeInputDeg : fThetaEpInputDeg;
    const G4double phiGenDeg = (pdg == 11) ? fPhiEeDeg : fPhiEpDeg;
    const G4double energyGenMeV = (pdg == 11) ? fEnergyEeMeV : fEnergyEpMeV;

    auto* analysis = G4AnalysisManager::Instance();

    // Tree 1: hits. One row per primary e-/e+ boundary crossing into a detector volume.
    // pdg = 11 is electron e-, pdg = -11 is positron e+.
    analysis->FillNtupleIColumn(1, 0, fEventID);
    analysis->FillNtupleIColumn(1, 1, fDataRow);
    analysis->FillNtupleIColumn(1, 2, track->GetTrackID());
    analysis->FillNtupleIColumn(1, 3, pdg);
    analysis->FillNtupleIColumn(1, 4, charge);
    analysis->FillNtupleIColumn(1, 5, detID);
    analysis->FillNtupleIColumn(1, 6, volumeID); // 0 = silicon, 1 = scintillator

    analysis->FillNtupleDColumn(1, 7, pos.x() / mm);
    analysis->FillNtupleDColumn(1, 8, pos.y() / mm);
    analysis->FillNtupleDColumn(1, 9, pos.z() / mm);
    analysis->FillNtupleDColumn(1, 10, pos.perp() / mm);

    analysis->FillNtupleDColumn(1, 11, pos.theta() / deg);
    analysis->FillNtupleDColumn(1, 12, phiHitDeg);
    analysis->FillNtupleDColumn(1, 13, pre->GetKineticEnergy() / MeV);
    analysis->FillNtupleDColumn(1, 14, pre->GetGlobalTime() / ns);
    analysis->FillNtupleDColumn(1, 15, mom.x());
    analysis->FillNtupleDColumn(1, 16, mom.y());
    analysis->FillNtupleDColumn(1, 17, mom.z());

    analysis->FillNtupleDColumn(1, 18, thetaGenDeg);
    analysis->FillNtupleDColumn(1, 19, phiGenDeg);
    analysis->FillNtupleDColumn(1, 20, energyGenMeV);

    analysis->AddNtupleRow(1);
}

void EventAction::EndOfEventAction(const G4Event*)
{
    auto* analysis = G4AnalysisManager::Instance();

    if (fDataRow < 0)
    {
        G4ExceptionDescription msg;
        msg << "No generated row was stored for event " << fEventID << ".\n"
            << "This would write -999 to the generated tree, so the run is stopped.\n"
            << "Rebuild the executable and make sure GeneratePrimaries() calls "
            << "EventAction::SetGeneratedKinematics().";
        G4Exception("EventAction::EndOfEventAction",
                    "X17GEN003",
                    FatalException,
                    msg);
    }

    // Tree 0: generated. One row per generated e-/e+ pair.
    analysis->FillNtupleIColumn(0, 0, fEventID);
    analysis->FillNtupleIColumn(0, 1, fDataRow);

    analysis->FillNtupleDColumn(0, 2, fThetaEEInputDeg);
    analysis->FillNtupleDColumn(0, 3, fThetaEeInputDeg);
    analysis->FillNtupleDColumn(0, 4, fEnergyEeMeV);
    analysis->FillNtupleDColumn(0, 5, fThetaEpInputDeg);
    analysis->FillNtupleDColumn(0, 6, fEnergyEpMeV);
    analysis->FillNtupleDColumn(0, 7, fPhi0Deg);
    analysis->FillNtupleDColumn(0, 8, fPhiEeDeg);
    analysis->FillNtupleDColumn(0, 9, fPhiEpDeg);

    analysis->FillNtupleDColumn(0, 10, fDirEe.x());
    analysis->FillNtupleDColumn(0, 11, fDirEe.y());
    analysis->FillNtupleDColumn(0, 12, fDirEe.z());
    analysis->FillNtupleDColumn(0, 13, fDirEp.x());
    analysis->FillNtupleDColumn(0, 14, fDirEp.y());
    analysis->FillNtupleDColumn(0, 15, fDirEp.z());
    analysis->FillNtupleDColumn(0, 16, OpeningAngleDeg(fDirEe, fDirEp));

    analysis->AddNtupleRow(0);

    Reset();
}
