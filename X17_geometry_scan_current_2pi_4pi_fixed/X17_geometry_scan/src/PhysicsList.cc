#include "PhysicsList.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
{
    SetDefaultCutValue(1.0 * mm);
}

void PhysicsList::ConstructParticle()
{
    // Particles needed by the generator.  No interaction processes are attached:
    // this is deliberately a transport-only setup, with no detector response effects.
    G4Gamma::GammaDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4OpticalPhoton::OpticalPhotonDefinition();
}

void PhysicsList::ConstructProcess()
{
    // Transportation only: straight propagation through the geometry.
    // No EM interactions, no multiple scattering, no energy loss, no annihilation.
    AddTransportation();
}

void PhysicsList::SetCuts()
{
    SetCutsWithDefault();
}
