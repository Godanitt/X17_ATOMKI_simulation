#include "EventAction.hh"
#include "X17EventInformation.hh"

#include "G4Event.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"

#include <algorithm>
#include <cmath>

EventAction::EventAction()
{
    ResetEvent();
}

void EventAction::ResetEvent()
{
    fEventID = -1;

    fThetaEEGenDeg = -999.0;
    fThetaEmGenDeg = -999.0;
    fThetaEpGenDeg = -999.0;
    fKinEmGenMeV = -999.0;
    fKinEpGenMeV = -999.0;

    fDirEmGen = G4ThreeVector(0.0, 0.0, 0.0);
    fDirEpGen = G4ThreeVector(0.0, 0.0, 0.0);

    fThetaEEIdealDeg = -999.0;
    fThetaEEIdealResDeg = -999.0;

    fHasEmSi = 0;
    fHasEpSi = 0;
    fHasEmScint = 0;
    fHasEpScint = 0;

    fDetEmSi = -1;
    fDetEpSi = -1;
    fDetEmScint = -1;
    fDetEpScint = -1;

    fPosEmSi = G4ThreeVector(0.0, 0.0, 0.0);
    fPosEpSi = G4ThreeVector(0.0, 0.0, 0.0);
    fPosEmScint = G4ThreeVector(0.0, 0.0, 0.0);
    fPosEpScint = G4ThreeVector(0.0, 0.0, 0.0);

    fThetaEEDetectorRawDeg = -999.0;

    fEdepEmScintMeV = 0.0;
    fEdepEpScintMeV = 0.0;
    fEdepAllScintMeV = 0.0;

    fNHitRows = 0;
    fNSiHitRows = 0;
    fNScintHitRows = 0;
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    ResetEvent();
    fEventID = event ? event->GetEventID() : -1;
}

void EventAction::SetGeneratedKinematics(G4double thetaEE_deg,
                                          G4double thetaEm_deg,
                                          G4double kinEm_MeV,
                                          G4double thetaEp_deg,
                                          G4double kinEp_MeV,
                                          const G4ThreeVector& dirEm,
                                          const G4ThreeVector& dirEp)
{
    fThetaEEGenDeg = thetaEE_deg;
    fThetaEmGenDeg = thetaEm_deg;
    fThetaEpGenDeg = thetaEp_deg;
    fKinEmGenMeV = kinEm_MeV;
    fKinEpGenMeV = kinEp_MeV;

    fDirEmGen = dirEm.unit();
    fDirEpGen = dirEp.unit();
}

G4double EventAction::OpeningAngleDeg(const G4ThreeVector& a,
                                      const G4ThreeVector& b) const
{
    if (a.mag2() <= 0.0 || b.mag2() <= 0.0) return -999.0;

    G4double cosang = a.unit().dot(b.unit());
    cosang = std::max(-1.0, std::min(1.0, cosang));

    return std::acos(cosang) / deg;
}

void EventAction::RecordDetectorStep(const G4Step* step,
                                     G4int detID,
                                     G4int volumeID,
                                     G4bool isBoundary)
{
    if (!step) return;

    const auto* track = step->GetTrack();
    const auto* pre = step->GetPreStepPoint();
    if (!track || !pre) return;

    const auto* particle = track->GetParticleDefinition();
    if (!particle) return;

    const G4int stepID = track->GetCurrentStepNumber();
    const G4int trackID = track->GetTrackID();
    const G4int parentID = track->GetParentID();
    const G4int pdg = particle->GetPDGEncoding();

    const G4double edepMeV = step->GetTotalEnergyDeposit() / MeV;
    const G4double stepLengthMm = step->GetStepLength() / mm;

    const G4ThreeVector pos = pre->GetPosition();

    const G4double xMm = pos.x() / mm;
    const G4double yMm = pos.y() / mm;
    const G4double zMm = pos.z() / mm;
    const G4double rMm = pos.mag() / mm;

    G4double thetaDeg = -999.0;
    G4double phiDeg = -999.0;

    if (rMm > 0.0)
    {
        thetaDeg = pos.theta() / deg;
        phiDeg = pos.phi() / deg;
        if (phiDeg < 0.0) phiDeg += 360.0;
    }

    const G4double kineticMeV = pre->GetKineticEnergy() / MeV;

    if (volumeID == 0) fNSiHitRows++;
    if (volumeID == 1) fNScintHitRows++;
    fNHitRows++;

    // Save the first boundary crossing as impact point. Prefer silicon for
    // angular reconstruction; scintillator positions are kept as a fallback.
    if (isBoundary && volumeID == 0)
    {
        if (pdg == 11 && fHasEmSi == 0)
        {
            fHasEmSi = 1;
            fDetEmSi = detID;
            fPosEmSi = pos;
        }
        else if (pdg == -11 && fHasEpSi == 0)
        {
            fHasEpSi = 1;
            fDetEpSi = detID;
            fPosEpSi = pos;
        }
    }

    if (isBoundary && volumeID == 1)
    {
        if (pdg == 11 && fHasEmScint == 0)
        {
            fHasEmScint = 1;
            fDetEmScint = detID;
            fPosEmScint = pos;
        }
        else if (pdg == -11 && fHasEpScint == 0)
        {
            fHasEpScint = 1;
            fDetEpScint = detID;
            fPosEpScint = pos;
        }
    }

    // The scintillator is the energy detector.
    if (volumeID == 1 && edepMeV > 0.0)
    {
        fEdepAllScintMeV += edepMeV;

        if (pdg == 11)       fEdepEmScintMeV += edepMeV;
        else if (pdg == -11) fEdepEpScintMeV += edepMeV;
    }

    auto* analysis = G4AnalysisManager::Instance();

    // Ntuple 1: hit/step-level information.
    analysis->FillNtupleIColumn(1, 0, fEventID);
    analysis->FillNtupleIColumn(1, 1, stepID);
    analysis->FillNtupleIColumn(1, 2, trackID);
    analysis->FillNtupleIColumn(1, 3, parentID);
    analysis->FillNtupleIColumn(1, 4, pdg);
    analysis->FillNtupleIColumn(1, 5, detID);
    analysis->FillNtupleIColumn(1, 6, volumeID);
    analysis->FillNtupleIColumn(1, 7, isBoundary ? 1 : 0);

    analysis->FillNtupleDColumn(1, 8, edepMeV);
    analysis->FillNtupleDColumn(1, 9, stepLengthMm);
    analysis->FillNtupleDColumn(1, 10, xMm);
    analysis->FillNtupleDColumn(1, 11, yMm);
    analysis->FillNtupleDColumn(1, 12, zMm);
    analysis->FillNtupleDColumn(1, 13, rMm);
    analysis->FillNtupleDColumn(1, 14, thetaDeg);
    analysis->FillNtupleDColumn(1, 15, phiDeg);
    analysis->FillNtupleDColumn(1, 16, kineticMeV);

    analysis->AddNtupleRow(1);
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Primary source of generator truth: G4Event user information.
    // This makes the ideal reconstruction independent of whether the generator
    // called SetGeneratedKinematics before or after BeginOfEventAction.
    if (event)
    {
        if (auto* info = dynamic_cast<X17EventInformation*>(event->GetUserInformation()))
        {
            fThetaEEGenDeg = info->thetaEE();
            fThetaEmGenDeg = info->thetaEm();
            fThetaEpGenDeg = info->thetaEp();
            fKinEmGenMeV = info->kinEm();
            fKinEpGenMeV = info->kinEp();
            fDirEmGen = info->dirEm().unit();
            fDirEpGen = info->dirEp().unit();
        }
    }

    fThetaEEIdealDeg = OpeningAngleDeg(fDirEmGen, fDirEpGen);
    fThetaEEIdealResDeg = fThetaEEIdealDeg - fThetaEEGenDeg;

    const G4bool hasBothSi = (fHasEmSi == 1 && fHasEpSi == 1);
    const G4bool hasBothScint = (fHasEmScint == 1 && fHasEpScint == 1);

    if (hasBothSi)
    {
        fThetaEEDetectorRawDeg = OpeningAngleDeg(fPosEmSi, fPosEpSi);
    }
    else if (hasBothScint)
    {
        fThetaEEDetectorRawDeg = OpeningAngleDeg(fPosEmScint, fPosEpScint);
    }

    const G4int pairGeomAccepted =
        ((hasBothSi || hasBothScint) ? 1 : 0);

    const G4int pairScintAccepted =
        (hasBothScint && fEdepEmScintMeV > 0.0 && fEdepEpScintMeV > 0.0) ? 1 : 0;

    auto* analysis = G4AnalysisManager::Instance();

    // Ntuple 0: event-level information.
    analysis->FillNtupleIColumn(0, 0, fEventID);

    analysis->FillNtupleDColumn(0, 1, fThetaEEGenDeg);
    analysis->FillNtupleDColumn(0, 2, fThetaEmGenDeg);
    analysis->FillNtupleDColumn(0, 3, fThetaEpGenDeg);
    analysis->FillNtupleDColumn(0, 4, fKinEmGenMeV);
    analysis->FillNtupleDColumn(0, 5, fKinEpGenMeV);

    analysis->FillNtupleDColumn(0, 6, fDirEmGen.x());
    analysis->FillNtupleDColumn(0, 7, fDirEmGen.y());
    analysis->FillNtupleDColumn(0, 8, fDirEmGen.z());
    analysis->FillNtupleDColumn(0, 9, fDirEpGen.x());
    analysis->FillNtupleDColumn(0, 10, fDirEpGen.y());
    analysis->FillNtupleDColumn(0, 11, fDirEpGen.z());

    analysis->FillNtupleDColumn(0, 12, fThetaEEIdealDeg);
    analysis->FillNtupleDColumn(0, 13, fThetaEEIdealResDeg);

    analysis->FillNtupleIColumn(0, 14, fHasEmSi);
    analysis->FillNtupleIColumn(0, 15, fHasEpSi);
    analysis->FillNtupleIColumn(0, 16, fHasEmScint);
    analysis->FillNtupleIColumn(0, 17, fHasEpScint);

    analysis->FillNtupleIColumn(0, 18, fDetEmSi);
    analysis->FillNtupleIColumn(0, 19, fDetEpSi);

    analysis->FillNtupleDColumn(0, 20, fPosEmSi.x() / mm);
    analysis->FillNtupleDColumn(0, 21, fPosEmSi.y() / mm);
    analysis->FillNtupleDColumn(0, 22, fPosEmSi.z() / mm);
    analysis->FillNtupleDColumn(0, 23, fPosEpSi.x() / mm);
    analysis->FillNtupleDColumn(0, 24, fPosEpSi.y() / mm);
    analysis->FillNtupleDColumn(0, 25, fPosEpSi.z() / mm);

    analysis->FillNtupleDColumn(0, 26, fThetaEEDetectorRawDeg);

    analysis->FillNtupleDColumn(0, 27, fEdepEmScintMeV);
    analysis->FillNtupleDColumn(0, 28, fEdepEpScintMeV);
    analysis->FillNtupleDColumn(0, 29, fEdepAllScintMeV);

    analysis->FillNtupleIColumn(0, 30, pairGeomAccepted);
    analysis->FillNtupleIColumn(0, 31, pairScintAccepted);

    analysis->FillNtupleIColumn(0, 32, fNHitRows);
    analysis->FillNtupleIColumn(0, 33, fNSiHitRows);
    analysis->FillNtupleIColumn(0, 34, fNScintHitRows);

    analysis->AddNtupleRow(0);
}
