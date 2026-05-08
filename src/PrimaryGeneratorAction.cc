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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction,
                                               const G4String& inputFile,
                                               const G4String& generationMode)
    : fGun(new G4ParticleGun(1)),
      fEventAction(eventAction)
{
    if (generationMode == "signal" || generationMode == "txt" || generationMode == "x17")
    {
        fMode = Mode::SignalFromTxt;
        LoadInputTable(inputFile);

        G4cout << "[generator] MODE = signal" << G4endl;
        G4cout << "[generator] INPUT TABLE = " << inputFile << G4endl;
        G4cout << "[generator] loaded " << fRows.size()
               << " rows from data/data_pair_creation.txt format:" << G4endl;
        G4cout << "[generator] col1 thetaEE_deg, col2 thetaEe_deg, col3 energyEe_MeV, "
               << "col4 thetaEp_deg, col5 energyEp_MeV" << G4endl;
    }
    else if (generationMode == "background" || generationMode == "bkg" || generationMode == "ipc")
    {
        fMode = Mode::SmoothBackground;

        G4cout << "[generator] MODE = background" << G4endl;
        G4cout << "[generator] smooth IPC-like toy background:" << G4endl;
        G4cout << "[generator]   total kinetic energy e-+e+ = "
               << fBackgroundTotalKineticMeV << " MeV" << G4endl;
        G4cout << "[generator]   per-particle kinetic minimum = "
               << fBackgroundMinKineticMeV << " MeV" << G4endl;
        G4cout << "[generator]   opening-angle scale = "
               << fBackgroundThetaScaleDeg << " deg" << G4endl;
        G4cout << "[generator] This is not a tuned IPC model; it is a smooth falling template "
               << "for detector/analysis studies." << G4endl;
    }
    else
    {
        G4ExceptionDescription msg;
        msg << "Unknown generation mode: " << generationMode << "\n"
            << "Allowed modes: signal, background.";
        G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction",
                    "X17GEN004",
                    FatalException,
                    msg);
    }
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
            << "The signal simulation intentionally has no fallback generator.\n"
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

G4ThreeVector PrimaryGeneratorAction::RandomIsotropicDirection() const
{
    const G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
    const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const G4double phi = CLHEP::twopi * G4UniformRand();

    return G4ThreeVector(sinTheta * std::cos(phi),
                         sinTheta * std::sin(phi),
                         cosTheta).unit();
}

G4ThreeVector PrimaryGeneratorAction::DirectionAtOpeningAngle(const G4ThreeVector& axis,
                                                               G4double openingAngle,
                                                               G4double rotationAngle) const
{
    const G4ThreeVector u = axis.unit();

    // Build a stable orthonormal basis around u.
    G4ThreeVector ref(0.0, 0.0, 1.0);
    if (std::abs(u.dot(ref)) > 0.95) ref = G4ThreeVector(1.0, 0.0, 0.0);

    const G4ThreeVector v = u.cross(ref).unit();
    const G4ThreeVector w = v.cross(u).unit();

    const G4ThreeVector transverse = (v * std::cos(rotationAngle) +
                                      w * std::sin(rotationAngle)).unit();
    return (u * std::cos(openingAngle) + transverse * std::sin(openingAngle)).unit();
}

G4double PrimaryGeneratorAction::SampleTruncatedExponentialDeg(G4double thetaMaxDeg,
                                                                G4double thetaScaleDeg) const
{
    // Sample theta in [0, thetaMax] with f(theta) ~ exp(-theta/thetaScale).
    // This is a toy IPC-like background: smooth and decreasing with opening angle.
    const G4double u = G4UniformRand();
    const G4double norm = 1.0 - std::exp(-thetaMaxDeg / thetaScaleDeg);
    return -thetaScaleDeg * std::log(1.0 - u * norm);
}

G4double PrimaryGeneratorAction::WrapPhiDeg(G4double phiDeg) const
{
    while (phiDeg < 0.0) phiDeg += 360.0;
    while (phiDeg >= 360.0) phiDeg -= 360.0;
    return phiDeg;
}

void PrimaryGeneratorAction::GenerateSignalFromTxt(G4Event* event, G4int eventID)
{
    // Use the supplied table deterministically: event 0 -> first data row,
    // event 1 -> second data row, and so on. If more events are requested than
    // rows exist, wrap around. No synthetic signal kinematic distribution is invented.
    const auto rowID = static_cast<std::size_t>(eventID) % fRows.size();
    const auto& row = fRows[rowID];

    // The table columns are:
    //   col1 thetaEE = opening angle between e- and e+
    //   col2 thetaEe = polar angle of the electron e-
    //   col3 energyEe = kinetic energy of the electron e-
    //   col4 thetaEp = polar angle of the positron e+
    //   col5 energyEp = kinetic energy of the positron e+
    //
    // We shoot the electron as e- and the positron as e+. The two particles are put
    // on opposite azimuthal sides of the same emission plane, with a random
    // global phi0 because the txt does not provide an absolute azimuth. This
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

void PrimaryGeneratorAction::GenerateSmoothBackground(G4Event* event, G4int eventID)
{
    // Smooth toy background:
    //   - opening angle is exponentially falling;
    //   - absolute orientation is isotropic;
    //   - total kinetic energy is fixed to the same scale as the txt sample;
    //   - energy sharing is random but keeps both particles above a small minimum.
    const G4double thetaEE_deg = SampleTruncatedExponentialDeg(180.0, fBackgroundThetaScaleDeg);
    const G4double thetaEE = thetaEE_deg * deg;

    const G4double available = fBackgroundTotalKineticMeV - 2.0 * fBackgroundMinKineticMeV;
    const G4double energyEe_MeV = fBackgroundMinKineticMeV + available * G4UniformRand();
    const G4double energyEp_MeV = fBackgroundTotalKineticMeV - energyEe_MeV;

    const G4ThreeVector dirEe = RandomIsotropicDirection();
    const G4double conePhi = CLHEP::twopi * G4UniformRand();
    const G4ThreeVector dirEp = DirectionAtOpeningAngle(dirEe, thetaEE, conePhi);

    const G4double thetaEe_deg = dirEe.theta() / deg;
    const G4double thetaEp_deg = dirEp.theta() / deg;
    const G4double phiEe_deg = WrapPhiDeg(dirEe.phi() / deg);
    const G4double phiEp_deg = WrapPhiDeg(dirEp.phi() / deg);
    const G4double phi0_deg = phiEe_deg;

    if (fEventAction)
    {
        fEventAction->SetGeneratedKinematics(eventID + 1,
                                             thetaEE_deg,
                                             thetaEe_deg,
                                             thetaEp_deg,
                                             energyEe_MeV,
                                             energyEp_MeV,
                                             phi0_deg,
                                             phiEe_deg,
                                             phiEp_deg,
                                             dirEe,
                                             dirEp);
    }

    const G4ThreeVector vertex(0.0, 0.0, 0.0);
    auto* table = G4ParticleTable::GetParticleTable();

    fGun->SetParticlePosition(vertex);

    fGun->SetParticleDefinition(table->FindParticle("e-"));
    fGun->SetParticleMomentumDirection(dirEe);
    fGun->SetParticleEnergy(energyEe_MeV * MeV);
    fGun->GeneratePrimaryVertex(event);

    fGun->SetParticleDefinition(table->FindParticle("e+"));
    fGun->SetParticleMomentumDirection(dirEp);
    fGun->SetParticleEnergy(energyEp_MeV * MeV);
    fGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const G4int eventID = event ? event->GetEventID() : 0;

    if (fMode == Mode::SignalFromTxt)
    {
        GenerateSignalFromTxt(event, eventID);
    }
    else
    {
        GenerateSmoothBackground(event, eventID);
    }
}
