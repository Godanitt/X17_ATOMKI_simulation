#include "PrimaryGeneratorAction.hh"

#include "EventAction.hh"

#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction,
                                               const G4String& inputFile)
    : fGun(new G4ParticleGun(1)),
      fEventAction(eventAction)
{
    LoadInputTable(inputFile);

    G4cout << "[generator] INPUT TABLE = " << inputFile << G4endl;
    G4cout << "[generator] loaded " << fRows.size()
           << " rows from data/data_pair_creation.txt format:" << G4endl;
    G4cout << "[generator] col1 thetaEE_deg, col2 thetaEe_deg, col3 energyEe_MeV, "
           << "col4 thetaEp_deg, col5 energyEp_MeV" << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

void PrimaryGeneratorAction::LoadInputTable(const G4String& inputFile)
{
    std::ifstream in(inputFile);
    if (!in.is_open())
    {
        G4ExceptionDescription msg;
        msg << "Cannot open input table: " << inputFile << "\n"
            << "The simulation intentionally has no fallback generator.\n"
            << "Run from the project root or pass -i /absolute/path/to/data_pair_creation.txt.";
        G4Exception("PrimaryGeneratorAction::LoadInputTable",
                    "X17GEN001",
                    FatalException,
                    msg);
    }

    std::string line;
    G4int lineNumber = 0;
    while (std::getline(in, line))
    {
        ++lineNumber;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        PairRow row;
        row.dataRow = lineNumber;

        // The five numbers are used directly from data/data_pair_creation.txt:
        //   thetaEE thetaEe energyEe thetaEp energyEp
        if (iss >> row.thetaEE_deg
                >> row.thetaEe_deg
                >> row.energyEe_MeV
                >> row.thetaEp_deg
                >> row.energyEp_MeV)
        {
            fRows.push_back(row);
        }
    }

    if (fRows.empty())
    {
        G4Exception("PrimaryGeneratorAction::LoadInputTable",
                    "X17GEN002",
                    FatalException,
                    "Input table is empty or badly formatted.");
    }
}

G4ThreeVector PrimaryGeneratorAction::DirectionFromPolarAzimuth(G4double theta,
                                                                  G4double phi) const
{
    return G4ThreeVector(std::sin(theta) * std::cos(phi),
                         std::sin(theta) * std::sin(phi),
                         std::cos(theta)).unit();
}

G4double PrimaryGeneratorAction::WrapPhiDeg(G4double phiDeg) const
{
    while (phiDeg < 0.0) phiDeg += 360.0;
    while (phiDeg >= 360.0) phiDeg -= 360.0;
    return phiDeg;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const G4int eventID = event ? event->GetEventID() : 0;

    // Use the supplied table deterministically: event 0 -> first data row,
    // event 1 -> second data row, and so on.  If more events are requested than
    // rows exist, wrap around.  No synthetic kinematic distribution is invented.
    const auto rowID = static_cast<std::size_t>(eventID) % fRows.size();
    const auto& row = fRows[rowID];

    // The table columns are:
    //   col1 thetaEE = opening angle between e- and e+
    //   col2 thetaEe = polar angle of the electron e-
    //   col3 energyEe = kinetic energy of the electron e-
    //   col4 thetaEp = polar angle of the positron e+
    //   col5 energyEp = kinetic energy of the positron e+
    //
    // We shoot the electron as e- and the positron as e+.  The two particles are put
    // on opposite azimuthal sides of the same emission plane, with a random
    // global phi0 because the txt does not provide an absolute azimuth.  This
    // preserves thetaEE = thetaEe + thetaEp for the supplied table.
    const G4double phi0 = CLHEP::twopi * G4UniformRand();
    const G4double phiEe = phi0;
    const G4double phiEp = phi0 + CLHEP::pi;

    const G4ThreeVector dirEe = DirectionFromPolarAzimuth(row.thetaEe_deg * deg, phiEe);
    const G4ThreeVector dirEp = DirectionFromPolarAzimuth(row.thetaEp_deg * deg, phiEp);

    if (fEventAction)
    {
        fEventAction->SetGeneratedKinematics(row.dataRow,
                                             row.thetaEE_deg,
                                             row.thetaEe_deg,
                                             row.thetaEp_deg,
                                             row.energyEe_MeV,
                                             row.energyEp_MeV,
                                             WrapPhiDeg(phi0 / deg),
                                             WrapPhiDeg(phiEe / deg),
                                             WrapPhiDeg(phiEp / deg),
                                             dirEe,
                                             dirEp);
    }

    const G4ThreeVector vertex(0.0, 0.0, 0.0);
    auto* table = G4ParticleTable::GetParticleTable();

    fGun->SetParticlePosition(vertex);

    fGun->SetParticleDefinition(table->FindParticle("e-"));
    fGun->SetParticleMomentumDirection(dirEe);
    fGun->SetParticleEnergy(row.energyEe_MeV * MeV);
    fGun->GeneratePrimaryVertex(event);

    fGun->SetParticleDefinition(table->FindParticle("e+"));
    fGun->SetParticleMomentumDirection(dirEp);
    fGun->SetParticleEnergy(row.energyEp_MeV * MeV);
    fGun->GeneratePrimaryVertex(event);
}
