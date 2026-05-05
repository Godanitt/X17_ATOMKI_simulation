#include "RunAction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"

RunAction::RunAction()
{
    auto* analysis = G4AnalysisManager::Instance();

    analysis->SetVerboseLevel(1);
    analysis->SetDefaultFileType("root");
    analysis->SetNtupleMerging(true);

    // ============================================================
    // NTUPLE 0: event-level truth, ideal reco, and compact detector summary
    // ============================================================
    analysis->CreateNtuple("events", "event-level generator and detector quantities");

    analysis->CreateNtupleIColumn("eventID");                    // 0

    analysis->CreateNtupleDColumn("thetaEE_gen_deg");            // 1
    analysis->CreateNtupleDColumn("thetaEm_gen_deg");            // 2
    analysis->CreateNtupleDColumn("thetaEp_gen_deg");            // 3
    analysis->CreateNtupleDColumn("kinEm_gen_MeV");              // 4
    analysis->CreateNtupleDColumn("kinEp_gen_MeV");              // 5

    analysis->CreateNtupleDColumn("dirEm_x");                    // 6
    analysis->CreateNtupleDColumn("dirEm_y");                    // 7
    analysis->CreateNtupleDColumn("dirEm_z");                    // 8
    analysis->CreateNtupleDColumn("dirEp_x");                    // 9
    analysis->CreateNtupleDColumn("dirEp_y");                    // 10
    analysis->CreateNtupleDColumn("dirEp_z");                    // 11

    analysis->CreateNtupleDColumn("thetaEE_ideal_deg");          // 12
    analysis->CreateNtupleDColumn("thetaEE_ideal_res_deg");      // 13

    analysis->CreateNtupleIColumn("hasEmSi");                    // 14
    analysis->CreateNtupleIColumn("hasEpSi");                    // 15
    analysis->CreateNtupleIColumn("hasEmScint");                 // 16
    analysis->CreateNtupleIColumn("hasEpScint");                 // 17

    analysis->CreateNtupleIColumn("detEmSi");                    // 18
    analysis->CreateNtupleIColumn("detEpSi");                    // 19

    analysis->CreateNtupleDColumn("xEmSi_mm");                   // 20
    analysis->CreateNtupleDColumn("yEmSi_mm");                   // 21
    analysis->CreateNtupleDColumn("zEmSi_mm");                   // 22
    analysis->CreateNtupleDColumn("xEpSi_mm");                   // 23
    analysis->CreateNtupleDColumn("yEpSi_mm");                   // 24
    analysis->CreateNtupleDColumn("zEpSi_mm");                   // 25

    analysis->CreateNtupleDColumn("thetaEE_detector_raw_deg");   // 26

    analysis->CreateNtupleDColumn("edepEmScint_MeV");            // 27
    analysis->CreateNtupleDColumn("edepEpScint_MeV");            // 28
    analysis->CreateNtupleDColumn("edepSumScint_MeV");           // 29

    analysis->CreateNtupleIColumn("pairGeomAccepted");           // 30
    analysis->CreateNtupleIColumn("pairScintAccepted");          // 31

    analysis->CreateNtupleIColumn("nHitRows");                   // 32
    analysis->CreateNtupleIColumn("nSiHitRows");                 // 33
    analysis->CreateNtupleIColumn("nScintHitRows");              // 34

    analysis->FinishNtuple();

    // ============================================================
    // NTUPLE 1: compact primary e-/e+ hit tree
    // ============================================================
    analysis->CreateNtuple("hits", "primary e-/e+ detector hit information");

    analysis->CreateNtupleIColumn("eventID");       // 0
    analysis->CreateNtupleIColumn("stepID");        // 1
    analysis->CreateNtupleIColumn("trackID");       // 2
    analysis->CreateNtupleIColumn("parentID");      // 3
    analysis->CreateNtupleIColumn("pdg");           // 4
    analysis->CreateNtupleIColumn("detID");         // 5
    analysis->CreateNtupleIColumn("volumeID");      // 6: 0 silicon, 1 scintillator
    analysis->CreateNtupleIColumn("isBoundary");    // 7

    analysis->CreateNtupleDColumn("edep_MeV");      // 8
    analysis->CreateNtupleDColumn("stepLength_mm"); // 9
    analysis->CreateNtupleDColumn("x_mm");          // 10
    analysis->CreateNtupleDColumn("y_mm");          // 11
    analysis->CreateNtupleDColumn("z_mm");          // 12
    analysis->CreateNtupleDColumn("r_mm");          // 13
    analysis->CreateNtupleDColumn("theta_deg");     // 14
    analysis->CreateNtupleDColumn("phi_deg");       // 15
    analysis->CreateNtupleDColumn("kinetic_MeV");   // 16

    analysis->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->OpenFile("x17_output.root");
    G4cout << "[RunAction] ROOT output file: x17_output.root" << G4endl;
}

void RunAction::EndOfRunAction(const G4Run*)
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->Write();
    analysis->CloseFile();
    G4cout << "[RunAction] ROOT file written." << G4endl;
}
