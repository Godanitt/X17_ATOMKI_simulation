#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
    : fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!fEventAction || !step) return;

    const auto* track = step->GetTrack();
    if (!track) return;

    const auto* particle = track->GetParticleDefinition();
    if (!particle) return;

    // For this exercise we only keep the primary e-/e+ detector response.
    const G4int pdg = particle->GetPDGEncoding();
    if (pdg != 11 && pdg != -11) return;
    if (track->GetParentID() != 0) return;

    const auto* pre = step->GetPreStepPoint();
    if (!pre) return;

    const auto* volume = pre->GetPhysicalVolume();
    if (!volume) return;

    const auto* logical = volume->GetLogicalVolume();
    if (!logical) return;

    const G4String lvName = logical->GetName();

    G4int volumeID = -1;
    if (lvName == "SiliconStripLV") volumeID = 0;
    if (lvName == "ScintillatorLV") volumeID = 1;

    if (volumeID < 0) return;

    const G4bool isBoundary = (pre->GetStepStatus() == fGeomBoundary);
    const G4double edep = step->GetTotalEnergyDeposit();

    // Keep boundary crossings for position reconstruction and non-zero
    // deposition steps for scintillator energy reconstruction.
    if (!isBoundary && edep <= 0.0) return;

    const auto touchable = pre->GetTouchableHandle();
    const G4int detID = touchable->GetCopyNumber();

    fEventAction->RecordDetectorStep(step, detID, volumeID, isBoundary);
}
