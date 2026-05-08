#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4ios.hh"

RunAction::RunAction(const G4String& outputFile)
    : fOutputFile(outputFile)
{
    auto* analysis = G4AnalysisManager::Instance();

    analysis->SetDefaultFileType("root");
    analysis->SetVerboseLevel(1);
    analysis->SetNtupleMerging(true);

    // Tree 0: generated. One row per event, using the row selected from
    // data/data_pair_creation.txt. The five *_txt_* values are direct copies
    // of the corresponding columns in that file:
    //   col1 thetaEE, col2 thetaEe, col3 energyEe, col4 thetaEp, col5 energyEp.
    analysis->CreateNtuple("generated", "generator-level e-/e+ pair kinematics from data/data_pair_creation.txt");
    analysis->CreateNtupleIColumn("eventID");
    analysis->CreateNtupleIColumn("dataRow");
    analysis->CreateNtupleDColumn("thetaEE_txt_deg");
    analysis->CreateNtupleDColumn("thetaEe_txt_deg");
    analysis->CreateNtupleDColumn("energyEe_txt_MeV");
    analysis->CreateNtupleDColumn("thetaEp_txt_deg");
    analysis->CreateNtupleDColumn("energyEp_txt_MeV");
    analysis->CreateNtupleDColumn("phi0_deg");
    analysis->CreateNtupleDColumn("phiEe_deg");
    analysis->CreateNtupleDColumn("phiEp_deg");
    analysis->CreateNtupleDColumn("dirEe_x");
    analysis->CreateNtupleDColumn("dirEe_y");
    analysis->CreateNtupleDColumn("dirEe_z");
    analysis->CreateNtupleDColumn("dirEp_x");
    analysis->CreateNtupleDColumn("dirEp_y");
    analysis->CreateNtupleDColumn("dirEp_z");
    analysis->CreateNtupleDColumn("thetaEE_generated_deg");
    analysis->FinishNtuple();

    // Tree 1: hits. One row per primary e-/e+ geometrical boundary crossing into
    // silicon/scintillator.  With the custom PhysicsList this is transport-only:
    // no smearing, no multiple scattering, no energy loss, no thresholds.
    // pdg = 11 is electron e-, pdg = -11 is positron e+.
    analysis->CreateNtuple("hits", "ideal transport-only detector boundary hits");
    analysis->CreateNtupleIColumn("eventID");
    analysis->CreateNtupleIColumn("dataRow");
    analysis->CreateNtupleIColumn("trackID");
    analysis->CreateNtupleIColumn("pdg");
    analysis->CreateNtupleIColumn("charge");
    analysis->CreateNtupleIColumn("detID");
    analysis->CreateNtupleIColumn("volumeID");       // 0 = SiliconStripLV, 1 = ScintillatorLV
    analysis->CreateNtupleDColumn("x_hit_mm");
    analysis->CreateNtupleDColumn("y_hit_mm");
    analysis->CreateNtupleDColumn("z_hit_mm");
    analysis->CreateNtupleDColumn("r_hit_mm");
    analysis->CreateNtupleDColumn("theta_hit_deg");
    analysis->CreateNtupleDColumn("phi_hit_deg");
    analysis->CreateNtupleDColumn("kinetic_hit_MeV");
    analysis->CreateNtupleDColumn("time_hit_ns");
    analysis->CreateNtupleDColumn("dir_hit_x");
    analysis->CreateNtupleDColumn("dir_hit_y");
    analysis->CreateNtupleDColumn("dir_hit_z");
    analysis->CreateNtupleDColumn("theta_gen_deg");
    analysis->CreateNtupleDColumn("phi_gen_deg");
    analysis->CreateNtupleDColumn("kinetic_gen_MeV");
    analysis->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->OpenFile(fOutputFile);
    G4cout << "[RunAction] ROOT output: " << fOutputFile << G4endl;
}

void RunAction::EndOfRunAction(const G4Run*)
{
    auto* analysis = G4AnalysisManager::Instance();
    analysis->Write();
    analysis->CloseFile();
    G4cout << "[RunAction] ROOT file written with exactly two trees: generated, hits." << G4endl;
}
