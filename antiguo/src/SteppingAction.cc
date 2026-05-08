#include "SteppingAction.hh"

#include "EventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

SteppingAction::SteppingAction(EventAction* eventAction)
    : fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    if (!fEventAction || !step) return;

    const auto* track = step->GetTrack();
    if (!track || track->GetParentID() != 0) return;

    const auto* particle = track->GetParticleDefinition();
    if (!particle) return;

    const auto pdg = particle->GetPDGEncoding();
    if (pdg != 11 && pdg != -11) return;

    const auto* pre = step->GetPreStepPoint();
    if (!pre || pre->GetStepStatus() != fGeomBoundary) return;

    const auto* volume = pre->GetPhysicalVolume();
    if (!volume || !volume->GetLogicalVolume()) return;

    const auto lvName = volume->GetLogicalVolume()->GetName();

    G4int volumeID = -1;
    if (lvName == "SiliconStripLV") volumeID = 0;
    else if (lvName == "ScintillatorLV") volumeID = 1;
    else return;

    const auto touchable = pre->GetTouchableHandle();
    const auto detID = touchable->GetCopyNumber();

    fEventAction->RecordHit(step, detID, volumeID);
}
