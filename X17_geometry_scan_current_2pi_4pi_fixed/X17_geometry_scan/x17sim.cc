#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4ios.hh"

#include <iostream>

namespace
{
void PrintUsage(const char* exe)
{
    std::cerr
        << "Usage:\n"
        << "  " << exe << " --vis [-i input.txt] [-o sampled.root] [--mode signal|background] [--geometry current|2pi|4pi|padplane]\n"
        << "  " << exe << " -m macro.mac [-i input.txt] [-o sampled.root] [--mode signal|background] [--geometry current|2pi|4pi|padplane]\n\n"
        << "Defaults:\n"
        << "  input    = data/data_pair_creation.txt\n"
        << "  output   = sampled.root\n"
        << "  mode     = signal\n"
        << "  geometry = current\n\n"
        << "Available geometries: " << DetectorConstruction::AvailableGeometries() << "\n";
}
}

int main(int argc, char** argv)
{
    G4String inputFile = "data/data_pair_creation.txt";
    G4String outputFile = "sampled.root";
    G4String macroFile;
    G4String generationMode = "signal";
    G4String geometryName = "current";
    G4bool visualMode = false;

    for (int i = 1; i < argc; ++i)
    {
        const G4String arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc)
        {
            inputFile = argv[++i];
        }
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc)
        {
            outputFile = argv[++i];
        }
        else if ((arg == "-m" || arg == "--macro") && i + 1 < argc)
        {
            macroFile = argv[++i];
        }
        else if ((arg == "--mode" || arg == "--sample") && i + 1 < argc)
        {
            generationMode = argv[++i];
        }
        else if ((arg == "--geometry" || arg == "--geom") && i + 1 < argc)
        {
            geometryName = argv[++i];
        }
        else if (arg == "--vis")
        {
            visualMode = true;
        }
        else if (arg == "-h" || arg == "--help")
        {
            PrintUsage(argv[0]);
            return 0;
        }
        else if (macroFile.empty())
        {
            macroFile = arg;
        }
        else
        {
            PrintUsage(argv[0]);
            return 1;
        }
    }

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    runManager->SetUserInitialization(new DetectorConstruction(geometryName));
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(inputFile, outputFile, generationMode));

    auto* visManager = new G4VisExecutive();
    visManager->Initialize();

    auto* uiManager = G4UImanager::GetUIpointer();

    G4cout << "[x17sim] Input file : " << inputFile << G4endl;
    G4cout << "[x17sim] Output file: " << outputFile << G4endl;
    G4cout << "[x17sim] Mode       : " << generationMode << G4endl;
    G4cout << "[x17sim] Geometry   : " << geometryName << G4endl;

    if (visualMode || macroFile.empty())
    {
        // Do not pass the program command-line options to G4UIExecutive.
        // Options such as --geometry or --vis are for this executable, not for
        // the Geant4 UI backend; passing them can prevent the viewer from opening
        // on some installations.
        char uiName[] = "x17sim";
        char* uiArgv[] = {uiName};
        int uiArgc = 1;
        auto* ui = new G4UIExecutive(uiArgc, uiArgv);
        uiManager->ApplyCommand("/control/execute macros/vis.mac");
        ui->SessionStart();
        delete ui;
    }
    else
    {
        uiManager->ApplyCommand("/control/execute " + macroFile);
    }

    delete visManager;
    delete runManager;

    return 0;
}
