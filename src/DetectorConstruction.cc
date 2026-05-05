#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <cmath>

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    // ----------------------------
    // Materials
    // ----------------------------
    auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    auto* targetMat = nist->FindOrBuildMaterial("G4_Li");
    auto* siliconMat = nist->FindOrBuildMaterial("G4_Si");
    auto* scintMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // ----------------------------
    // World
    // ----------------------------
    auto* solidWorld = new G4Box("World", 70.0 * cm, 70.0 * cm, 70.0 * cm);

    auto* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "WorldLV");

    auto* physWorld = new G4PVPlacement(nullptr,
                                        G4ThreeVector(),
                                        logicWorld,
                                        "WorldPV",
                                        nullptr,
                                        false,
                                        0,
                                        true);

    auto* worldVis = new G4VisAttributes();
    worldVis->SetVisibility(false);
    logicWorld->SetVisAttributes(worldVis);

    // ----------------------------
    // Thin lithium target. Beam axis = Z.
    // ----------------------------
    auto* solidTarget = new G4Tubs("LiTarget",
                                   0.0,
                                   5.0 * mm,
                                   0.05 * mm,
                                   0.0,
                                   360.0 * deg);

    auto* logicTarget = new G4LogicalVolume(solidTarget, targetMat, "LiTargetLV");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0.0, 0.0, 0.0),
                      logicTarget,
                      "LiTargetPV",
                      logicWorld,
                      false,
                      0,
                      true);

    auto* targetVis = new G4VisAttributes(G4Colour(0.9, 0.2, 0.2));
    targetVis->SetForceSolid(true);
    logicTarget->SetVisAttributes(targetVis);

    // ----------------------------
    // ATOMKI-like hexagonal telescope array
    // ----------------------------
    // Each telescope has:
    //   1) an inner silicon strip plane for the impact point;
    //   2) an outer plastic scintillator for the deposited energy.
    //
    // Local axes of each telescope:
    //   local x = radial thickness
    //   local y = tangential width
    //   local z = height along the beam axis
    // ----------------------------
    const G4int nDet = 6;

    const G4double scintHalfThickness = 1.0 * cm;
    const G4double scintHalfWidth     = 4.1 * cm;
    const G4double scintHalfHeight    = 4.3 * cm;

    const G4double siHalfThickness = 0.25 * mm;  // 0.5 mm thick detector
    const G4double siHalfWidth     = 4.1 * cm;
    const G4double siHalfHeight    = 4.3 * cm;

    const G4double gapSiScint = 2.0 * mm;

    const G4double hexSide = 2.0 * scintHalfWidth;
    const G4double hexApothem = 0.5 * std::sqrt(3.0) * hexSide;

    // Radius of scintillator centers. The small factor opens a visible gap.
    const G4double scintRadius = 1.03 * (hexApothem + scintHalfThickness);
    const G4double siliconRadius =
        scintRadius - scintHalfThickness - gapSiScint - siHalfThickness;

    auto* solidScint = new G4Box("Scintillator",
                                  scintHalfThickness,
                                  scintHalfWidth,
                                  scintHalfHeight);

    auto* logicScint = new G4LogicalVolume(solidScint,
                                           scintMat,
                                           "ScintillatorLV");

    auto* solidSilicon = new G4Box("SiliconStrip",
                                   siHalfThickness,
                                   siHalfWidth,
                                   siHalfHeight);

    auto* logicSilicon = new G4LogicalVolume(solidSilicon,
                                             siliconMat,
                                             "SiliconStripLV");

    auto* scintVis = new G4VisAttributes(G4Colour(0.1, 0.4, 1.0, 0.28));
    scintVis->SetForceSolid(true);
    logicScint->SetVisAttributes(scintVis);

    auto* siliconVis = new G4VisAttributes(G4Colour(0.55, 0.15, 1.0, 0.75));
    siliconVis->SetForceSolid(true);
    logicSilicon->SetVisAttributes(siliconVis);

    for (G4int i = 0; i < nDet; ++i)
    {
        const G4double phi = i * 360.0 * deg / nDet;

        auto* rot = new G4RotationMatrix();
        rot->rotateZ(phi);

        const G4double xs = scintRadius * std::cos(phi);
        const G4double ys = scintRadius * std::sin(phi);

        const G4double xi = siliconRadius * std::cos(phi);
        const G4double yi = siliconRadius * std::sin(phi);

        new G4PVPlacement(rot,
                          G4ThreeVector(xi, yi, 0.0),
                          logicSilicon,
                          "SiliconStripPV",
                          logicWorld,
                          false,
                          i,
                          true);

        new G4PVPlacement(rot,
                          G4ThreeVector(xs, ys, 0.0),
                          logicScint,
                          "ScintillatorPV",
                          logicWorld,
                          false,
                          i,
                          true);
    }

    return physWorld;
}
