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
#include <cstdlib>

namespace
{
    G4double Clamp(G4double x, G4double lo, G4double hi)
    {
        return std::max(lo, std::min(hi, x));
    }

    G4double ReadEnvDouble(const char* name, G4double fallback)
    {
        const char* value = std::getenv(name);
        if (!value) return fallback;

        char* end = nullptr;
        const double parsed = std::strtod(value, &end);
        if (end == value || !std::isfinite(parsed)) return fallback;

        return parsed;
    }
}

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction,
                                               const G4String& filename)
    : fGun(nullptr),
      fEventAction(eventAction),
      fSignalFraction(0.10),
      fIPCThetaMinDeg(20.0),
      fIPCThetaMaxDeg(180.0),
      fIPCThetaSlopePerDeg(0.018),
      fIPCEnergyAsymmetryMax(0.85),
      fIPCEnergyScaleSigma(0.04)
{
    fGun = new G4ParticleGun(1);
    LoadTable(filename);
    ConfigureFromEnvironment();

    G4cout << "[X17 generator] loaded " << fTable.size()
           << " X17-like rows from " << filename << G4endl;

    G4cout << "[X17 generator] event mixture: IPC-like smooth background = "
           << 100.0 * (1.0 - fSignalFraction) << " %, X17-like signal = "
           << 100.0 * fSignalFraction << " %" << G4endl;

    G4cout << "[X17 generator] IPC model: theta_ee in ["
           << fIPCThetaMinDeg << ", " << fIPCThetaMaxDeg
           << "] deg with dN/dtheta proportional to exp(-"
           << fIPCThetaSlopePerDeg << " * theta_ee/deg)" << G4endl;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

void PrimaryGeneratorAction::ConfigureFromEnvironment()
{
    // Runtime knobs without adding extra Geant4 messenger code. Example:
    //   X17_SIGNAL_FRACTION=0.05 X17_IPC_THETA_SLOPE_PER_DEG=0.02 ./build/x17sim logs/run_default.mac
    fSignalFraction = Clamp(ReadEnvDouble("X17_SIGNAL_FRACTION", fSignalFraction), 0.0, 1.0);

    fIPCThetaMinDeg = Clamp(ReadEnvDouble("X17_IPC_THETA_MIN_DEG", fIPCThetaMinDeg), 0.0, 180.0);
    fIPCThetaMaxDeg = Clamp(ReadEnvDouble("X17_IPC_THETA_MAX_DEG", fIPCThetaMaxDeg), 0.0, 180.0);
    if (fIPCThetaMaxDeg <= fIPCThetaMinDeg)
    {
        fIPCThetaMinDeg = 20.0;
        fIPCThetaMaxDeg = 180.0;
    }

    fIPCThetaSlopePerDeg = std::max(0.0, ReadEnvDouble("X17_IPC_THETA_SLOPE_PER_DEG",
                                                       fIPCThetaSlopePerDeg));
    fIPCEnergyAsymmetryMax = Clamp(ReadEnvDouble("X17_IPC_ENERGY_ASYMMETRY_MAX",
                                                 fIPCEnergyAsymmetryMax),
                                   0.0, 0.98);
    fIPCEnergyScaleSigma = Clamp(ReadEnvDouble("X17_IPC_ENERGY_SCALE_SIGMA",
                                               fIPCEnergyScaleSigma),
                                 0.0, 0.30);
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

PrimaryGeneratorAction::PairEvent PrimaryGeneratorAction::SampleSignalPair() const
{
    const G4int irow = G4RandFlat::shootInt(fTable.size());
    return fTable[irow];
}

G4double PrimaryGeneratorAction::SampleTruncatedExponentialDeg(G4double thetaMinDeg,
                                                               G4double thetaMaxDeg,
                                                               G4double slopePerDeg) const
{
    const G4double width = thetaMaxDeg - thetaMinDeg;
    if (width <= 0.0) return thetaMinDeg;

    const G4double u = G4UniformRand();

    if (slopePerDeg <= 0.0)
    {
        return thetaMinDeg + width * u;
    }

    const G4double norm = 1.0 - std::exp(-slopePerDeg * width);
    return thetaMinDeg - std::log(1.0 - u * norm) / slopePerDeg;
}

PrimaryGeneratorAction::PairEvent PrimaryGeneratorAction::SampleIPCBackgroundPair() const
{
    // IPC-like background model:
    //   * no resonant bump;
    //   * aperture angle is smooth and decreasing;
    //   * total kinetic-energy scale is borrowed from the supplied pair table so
    //     the detector response stays in the same dynamic range as the signal;
    //   * lepton energy sharing is broad, as expected for continuum pairs.
    PairEvent ev = SampleSignalPair();

    ev.thetaEE_deg = SampleTruncatedExponentialDeg(fIPCThetaMinDeg,
                                                   fIPCThetaMaxDeg,
                                                   fIPCThetaSlopePerDeg);

    const G4double templateTotal = std::max(0.10, ev.kineticElectron_MeV + ev.kineticPositron_MeV);
    const G4double scale = std::max(0.20, G4RandGauss::shoot(1.0, fIPCEnergyScaleSigma));
    const G4double total = templateTotal * scale;

    const G4double y = G4RandFlat::shoot(-fIPCEnergyAsymmetryMax,
                                         +fIPCEnergyAsymmetryMax);

    ev.kineticElectron_MeV = std::max(0.02, 0.5 * total * (1.0 + y));
    ev.kineticPositron_MeV = std::max(0.02, 0.5 * total * (1.0 - y));

    // These fields are not used for the transverse-plane generator below, but
    // keeping physically plausible values makes the truth tree self-describing.
    ev.thetaElectron_deg = 90.0;
    ev.thetaPositron_deg = 90.0;

    return ev;
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
    const G4bool makeSignal = (G4UniformRand() < fSignalFraction);
    const G4int componentID = makeSignal ? kX17Signal : kIPCBackground;
    const PairEvent row = makeSignal ? SampleSignalPair() : SampleIPCBackgroundPair();

    // ------------------------------------------------------------
    // Reconstruction-first generator.
    // ------------------------------------------------------------
    // Both leptons are generated in the transverse plane, where the
    // ATOMKI-like telescope array sits. The X17 component uses the supplied
    // table value thetaEE_deg. The IPC component replaces that value with a
    // smooth, monotonically decreasing continuum distribution.
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
                                             dirPositron,
                                             componentID);
    }

    event->SetUserInformation(
        new X17EventInformation(row.thetaEE_deg,
                                thetaLab / deg,
                                row.kineticElectron_MeV,
                                thetaLab / deg,
                                row.kineticPositron_MeV,
                                dirElectron,
                                dirPositron,
                                componentID));

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
               << "component = " << (componentID == kX17Signal ? "X17-like signal" : "IPC-like background") << ", "
               << "theta_ee(input) = " << row.thetaEE_deg << " deg, "
               << "theta_ee(generated) = " << opening << " deg, "
               << "theta_lab(e-) = " << thetaLab / deg << " deg, "
               << "theta_lab(e+) = " << thetaLab / deg << " deg, "
               << "T(e-) = " << row.kineticElectron_MeV << " MeV, "
               << "T(e+) = " << row.kineticPositron_MeV << " MeV"
               << G4endl;
    }
}
