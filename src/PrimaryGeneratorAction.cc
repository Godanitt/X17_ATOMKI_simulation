#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "X17EventInformation.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4Exception.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction,
                                               const G4String& filename)
    : fGun(nullptr),
      fEventAction(eventAction)
{
    fGun = new G4ParticleGun(1);
    LoadTable(filename);

    G4cout << "[X17 generator] loaded " << fTable.size()
           << " rows from " << filename << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

void PrimaryGeneratorAction::LoadTable(const G4String& filename)
{
    std::ifstream input(filename);

    if (!input.is_open())
    {
        G4ExceptionDescription msg;
        msg << "Cannot open " << filename << ".\n"
            << "Run from the project root or make sure data/data_pair_creation.txt exists.";

        G4Exception("PrimaryGeneratorAction::LoadTable",
                    "X17GEN001",
                    FatalException,
                    msg);
    }

    std::string line;

    while (std::getline(input, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        PairEvent ev;

        if (iss >> ev.thetaEE_deg
                >> ev.thetaElectron_deg
                >> ev.kineticElectron_MeV
                >> ev.thetaPositron_deg
                >> ev.kineticPositron_MeV)
        {
            fTable.push_back(ev);
        }
    }

    if (fTable.empty())
    {
        G4Exception("PrimaryGeneratorAction::LoadTable",
                    "X17GEN002",
                    FatalException,
                    "The input table is empty or badly formatted.");
    }
}

G4ThreeVector PrimaryGeneratorAction::DirectionFromThetaPhi(G4double theta,
                                                            G4double phi) const
{
    return G4ThreeVector(std::sin(theta) * std::cos(phi),
                         std::sin(theta) * std::sin(phi),
                         std::cos(theta));
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const G4int irow = G4RandFlat::shootInt(fTable.size());
    const PairEvent& row = fTable[irow];

    // ------------------------------------------------------------
    // Reconstruction-first generator.
    // ------------------------------------------------------------
    // For this stage we want a clean angular reconstruction before
    // introducing detector effects. We therefore generate both leptons
    // in the transverse plane, where the ATOMKI-like telescope array sits,
    // and use the table value thetaEE_deg as their true 3D opening angle.
    //
    // Detector effects are studied later in the offline analysis by using
    // silicon hit positions and scintillator energy deposits.
    // ------------------------------------------------------------
    const G4double thetaLab = 90.0 * deg;
    const G4double thetaEE = row.thetaEE_deg * deg;

    const G4double phiPair = CLHEP::twopi * G4UniformRand();
    const G4double phiElectron = phiPair - 0.5 * thetaEE;
    const G4double phiPositron = phiPair + 0.5 * thetaEE;

    const G4ThreeVector dirElectron =
        DirectionFromThetaPhi(thetaLab, phiElectron).unit();

    const G4ThreeVector dirPositron =
        DirectionFromThetaPhi(thetaLab, phiPositron).unit();

    const G4double kineticElectron = row.kineticElectron_MeV * MeV;
    const G4double kineticPositron = row.kineticPositron_MeV * MeV;

    // Store truth twice for robustness:
    //   1) directly in the EventAction used by this worker,
    //   2) inside the G4Event as user information, read back at EndOfEventAction.
    // This protects the ROOT truth branches against action-order or MT surprises.
    if (fEventAction)
    {
        fEventAction->SetGeneratedKinematics(row.thetaEE_deg,
                                             thetaLab / deg,
                                             row.kineticElectron_MeV,
                                             thetaLab / deg,
                                             row.kineticPositron_MeV,
                                             dirElectron,
                                             dirPositron);
    }

    event->SetUserInformation(
        new X17EventInformation(row.thetaEE_deg,
                                thetaLab / deg,
                                row.kineticElectron_MeV,
                                thetaLab / deg,
                                row.kineticPositron_MeV,
                                dirElectron,
                                dirPositron));

    const G4ThreeVector vertex(0.0, 0.0, 0.0);
    auto* particleTable = G4ParticleTable::GetParticleTable();

    fGun->SetParticleDefinition(particleTable->FindParticle("e-"));
    fGun->SetParticlePosition(vertex);
    fGun->SetParticleMomentumDirection(dirElectron);
    fGun->SetParticleEnergy(kineticElectron);
    fGun->GeneratePrimaryVertex(event);

    fGun->SetParticleDefinition(particleTable->FindParticle("e+"));
    fGun->SetParticlePosition(vertex);
    fGun->SetParticleMomentumDirection(dirPositron);
    fGun->SetParticleEnergy(kineticPositron);
    fGun->GeneratePrimaryVertex(event);

    if (event->GetEventID() < 5)
    {
        G4double cosOpening = dirElectron.dot(dirPositron);
        cosOpening = std::max(-1.0, std::min(1.0, cosOpening));
        const G4double opening = std::acos(cosOpening) / deg;

        G4cout << "[event " << event->GetEventID() << "] "
               << "theta_ee(table) = " << row.thetaEE_deg << " deg, "
               << "theta_ee(generated) = " << opening << " deg, "
               << "theta_lab(e-) = " << thetaLab / deg << " deg, "
               << "theta_lab(e+) = " << thetaLab / deg << " deg, "
               << "T(e-) = " << row.kineticElectron_MeV << " MeV, "
               << "T(e+) = " << row.kineticPositron_MeV << " MeV"
               << G4endl;
    }
}
