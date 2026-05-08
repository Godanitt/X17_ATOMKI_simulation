#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <string>
#include <vector>

namespace
{
struct DetectorLogicals
{
    G4LogicalVolume* world = nullptr;
    G4LogicalVolume* silicon = nullptr;
    G4LogicalVolume* scintillator = nullptr;

    G4double scintHalfThickness = 1.0 * cm;
    G4double scintHalfWidth = 4.1 * cm;
    G4double scintHalfHeight = 4.3 * cm;

    G4double siHalfThickness = 0.25 * mm;
    G4double siHalfWidth = 4.1 * cm;
    G4double siHalfHeight = 4.3 * cm;

    G4double gapSiScint = 2.0 * mm;
};

std::string lower_copy(const G4String& in)
{
    std::string s = in;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

void PlaceTelescopePair(const DetectorLogicals& det,
                        const G4ThreeVector& siliconPosition,
                        const G4ThreeVector& scintillatorPosition,
                        G4RotationMatrix* rotation,
                        G4int copyNo,
                        const G4String& namePrefix)
{
    new G4PVPlacement(rotation,
                      siliconPosition,
                      det.silicon,
                      namePrefix + "SiliconPV",
                      det.world,
                      false,
                      copyNo,
                      true);

    new G4PVPlacement(rotation,
                      scintillatorPosition,
                      det.scintillator,
                      namePrefix + "ScintillatorPV",
                      det.world,
                      false,
                      copyNo,
                      true);
}


G4RotationMatrix* RotationWithLocalXAlong(const G4ThreeVector& direction)
{
    // Detector boxes are defined with their thin dimension along local X.
    // Therefore the local X axis is the normal to the detector plane.
    //
    // The correct telescope orientation is:
    //   local +X  -> radial outward direction from the target
    //   local Y,Z -> tangent directions to the sphere/cylinder
    //
    // Using explicit axes is safer than composing rotateY/rotateZ calls, because
    // the latter is easy to get wrong and produced visibly tilted modules in the
    // 2pi/4pi views.
    const G4ThreeVector xAxis = direction.unit();

    // Reference axis used to define a stable tangent direction.  If the detector
    // points close to the beam axis, use global Y instead of global Z to avoid a
    // nearly-zero cross product.
    const G4ThreeVector zRef(0.0, 0.0, 1.0);
    const G4ThreeVector yRef(0.0, 1.0, 0.0);
    const G4ThreeVector ref = (std::abs(xAxis.dot(zRef)) < 0.95) ? zRef : yRef;

    const G4ThreeVector yAxis = ref.cross(xAxis).unit();
    const G4ThreeVector zAxis = xAxis.cross(yAxis).unit();

    auto* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);
    return rot;
}

void PlaceDirectionalTelescope(const DetectorLogicals& det,
                               const G4ThreeVector& direction,
                               G4double siliconRadius,
                               G4int copyNo,
                               const G4String& namePrefix)
{
    const G4ThreeVector u = direction.unit();
    const G4double scintRadius =
        siliconRadius + det.siHalfThickness + det.gapSiScint + det.scintHalfThickness;

    PlaceTelescopePair(det,
                       siliconRadius * u,
                       scintRadius * u,
                       RotationWithLocalXAlong(u),
                       copyNo,
                       namePrefix);
}

G4int BuildHemisphereShell(const DetectorLogicals& det,
                           G4int firstCopyNo,
                           G4double siliconRadius,
                           G4bool downstreamHemisphere,
                           const G4String& namePrefix,
                           G4double rotXDeg = 0.0,
                           G4double rotYDeg = 0.0,
                           G4double rotZDeg = 0.0)
{
    struct RingSpec
    {
        G4double thetaFromAxisDeg;
        G4int nPhi;
        G4double phiOffsetDeg;
    };

    const std::vector<RingSpec> rings = {
        {28.0, 6,  0.0},
        {55.0, 10, 18.0},
        {90.0, 14,  0.0}
    };

    G4int copyNo = firstCopyNo;

    for (const auto& ring : rings)
    {
        const G4double theta = ring.thetaFromAxisDeg * deg;

        for (G4int i = 0; i < ring.nPhi; ++i)
        {

            G4RotationMatrix globalRot;
            globalRot.rotateX(rotXDeg * deg);
            globalRot.rotateY(rotYDeg * deg);
            globalRot.rotateZ(rotZDeg * deg);

            const G4double phi = (ring.phiOffsetDeg + i * 360.0 / ring.nPhi) * deg;

            const G4ThreeVector direction0(std::sin(theta) * std::cos(phi),
                                           std::sin(theta) * std::sin(phi),
                                           std::cos(theta));

            const G4ThreeVector direction = (globalRot * direction0).unit();

            PlaceDirectionalTelescope(det,
                                      direction,
                                      siliconRadius,
                                      copyNo,
                                      namePrefix);

            ++copyNo;
        }
    }

    return copyNo;
}

void BuildCurrentAtomkiLikeGeometry(const DetectorLogicals& det)
{
    // Original simplified ATOMKI-like hexagonal telescope array.
    const G4int nDet = 6;

    const G4double hexSide = 2.0 * det.scintHalfWidth;
    const G4double hexApothem = 0.5 * std::sqrt(3.0) * hexSide;

    const G4double scintRadius = 1.03 * (hexApothem + det.scintHalfThickness);
    const G4double siliconRadius =
        scintRadius - det.scintHalfThickness - det.gapSiScint - det.siHalfThickness;

    for (G4int i = 0; i < nDet; ++i)
    {
        const G4double phi = i * 360.0 * deg / nDet;

        const G4double xs = 1.1 * scintRadius * std::cos(phi);
        const G4double ys = 1.1 * scintRadius * std::sin(phi);

        const G4double xi = 1.1 * siliconRadius * std::cos(phi);
        const G4double yi = 1.1 * siliconRadius * std::sin(phi);

        const G4ThreeVector u(std::cos(phi), std::sin(phi), 0.0);

        PlaceTelescopePair(det,
                           G4ThreeVector(xi, yi, 0.0),
                           G4ThreeVector(xs, ys, 0.0),
                           RotationWithLocalXAlong(u),
                           i,
                           "Current");
    }
}

void BuildBarrelRing(const DetectorLogicals& det,
                     G4int nDetectors,
                     G4double siliconRadius,
                     G4int firstCopyNo,
                     const G4String& namePrefix,
                     G4double zCenter = 0.0)
{
    const G4double scintRadius =
        siliconRadius + det.siHalfThickness + det.gapSiScint + det.scintHalfThickness;

    for (G4int i = 0; i < nDetectors; ++i)
    {
        const G4double phi = i * 360.0 * deg / nDetectors;

        const G4ThreeVector u(std::cos(phi), std::sin(phi), 0.0);
        const G4ThreeVector siliconPosition = siliconRadius * u + G4ThreeVector(0.0, 0.0, zCenter);
        const G4ThreeVector scintPosition = scintRadius * u + G4ThreeVector(0.0, 0.0, zCenter);

        PlaceTelescopePair(det,
                           siliconPosition,
                           scintPosition,
                           RotationWithLocalXAlong(u),
                           firstCopyNo + i,
                           namePrefix);
    }
}

G4int BuildPadPlane(const DetectorLogicals& det,
                    G4int firstCopyNo,
                    G4double zSilicon,
                    G4int nx,
                    G4int ny,
                    G4double pitchX,
                    G4double pitchY,
                    const G4String& namePrefix)
{
    const G4bool downstream = zSilicon >= 0.0;
    const G4double sign = downstream ? 1.0 : -1.0;

    const G4double zScint = zSilicon + sign * (det.siHalfThickness + det.gapSiScint + det.scintHalfThickness);

    // Ángulo de giro alrededor del eje Z
    const G4double alpha = 90.0 * deg;
    const G4double ca = std::cos(alpha);
    const G4double sa = std::sin(alpha);

    G4int copyNo = firstCopyNo;
    for (G4int ix = 0; ix < nx; ++ix)
    {
        for (G4int iy = 0; iy < ny; ++iy)
        {
            const G4double x0 = (ix - 0.5 * (nx - 1)) * pitchX;
            const G4double y0 = (iy - 0.5 * (ny - 1)) * pitchY;

            // Hueco central para el haz
            if (std::abs(x0) < 1.e-9 * mm && std::abs(y0) < 1.e-9 * mm)
            {
                continue;
            }

            // Giro alrededor del eje Z
            const G4double x = x0 * ca - y0 * sa;
            const G4double y = x0 * sa + y0 * ca;

            const G4ThreeVector normal(0.0, 0.0, sign);

            auto* rot = RotationWithLocalXAlong(normal);

            // Esto solo hace falta si el detector no es cuadrado
            // o si quieres que sus ejes locales también roten en el plano.
            rot->rotateZ(alpha);

            PlaceTelescopePair(det,
                               G4ThreeVector(x, y, zSilicon),
                               G4ThreeVector(x, y, zScint),
                               rot,
                               copyNo,
                               namePrefix);
            ++copyNo;
        }
    }

    return copyNo;
}

void BuildTwoPiCoverageGeometry(const DetectorLogicals& det)
{
    // Semispherical 2pi coverage: detector telescopes distributed on the downstream
    // hemisphere (+Z side of the target). This is intentionally predefined: no
    // radius or multiplicity has to be provided from the macro or command line.
    BuildHemisphereShell(det,
                         0,
                         22.0 * cm,
                         true,
                         "TwoPiHemisphere");
}

void BuildFourPiCoverageGeometry(const DetectorLogicals& det)
{
    // Full 4pi-inspired layout: two mirrored semispherical shells around the target.
    // It reuses the same module technology and the same sensitive volumes as the
    // current and 2pi configurations.
    G4int nextCopyNo = BuildHemisphereShell(det,
                                            0,
                                            22.0 * cm,
                                            true,
                                            "FourPiForwardHemisphere");

    BuildHemisphereShell(det,
                         nextCopyNo,
                         22.0 * cm,
                         false,
                         "FourPiBackwardHemisphere");
}

void BuildForwardPadPlaneGeometry(const DetectorLogicals& det)
{
    // Plane of pads downstream of the target, perpendicular to the beam axis (+z).
    BuildPadPlane(det,
                  0,
                  +18.0 * cm,
                  5,
                  5,
                  9.0 * cm,
                  9.0 * cm,
                  "ForwardPadPlane");
}
}

DetectorConstruction::DetectorConstruction(const G4String& geometryName)
{
    DefineCommands();
    SetGeometry(geometryName);
}

DetectorConstruction::~DetectorConstruction()
{
    delete fMessenger;
}

G4String DetectorConstruction::AvailableGeometries()
{
    return "current, 2pi, 4pi, padplane";
}

void DetectorConstruction::DefineCommands()
{
    fMessenger = new G4GenericMessenger(this,
                                        "/x17/geometry/",
                                        "X17 detector geometry control");

    fMessenger->DeclareMethod("select",
                              &DetectorConstruction::SetGeometry,
                              "Select one predefined geometry: current, 2pi, 4pi, padplane. Must be called before /run/initialize.");
}

void DetectorConstruction::SetGeometry(const G4String& geometryName)
{
    const auto key = lower_copy(geometryName);

    if (key == "current" || key == "atomki" || key == "baseline")
    {
        fGeometry = GeometryType::Current;
        fGeometryName = "current";
    }
    else if (key == "2pi" || key == "coverage2pi" || key == "two_pi")
    {
        fGeometry = GeometryType::Coverage2Pi;
        fGeometryName = "2pi";
    }
    else if (key == "4pi" || key == "coverage4pi" || key == "four_pi")
    {
        fGeometry = GeometryType::Coverage4Pi;
        fGeometryName = "4pi";
    }
    else if (key == "padplane" || key == "pad_plane" || key == "pads" || key == "frontpads")
    {
        fGeometry = GeometryType::PadPlane;
        fGeometryName = "padplane";
    }
    else
    {
        G4ExceptionDescription msg;
        msg << "Unknown geometry '" << geometryName << "'. Available geometries: "
            << AvailableGeometries();
        G4Exception("DetectorConstruction::SetGeometry",
                    "X17GEOM001",
                    FatalException,
                    msg);
    }

    G4cout << "[DetectorConstruction] Selected geometry: " << fGeometryName << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
    auto* targetMat = nist->FindOrBuildMaterial("G4_Li");
    auto* siliconMat = nist->FindOrBuildMaterial("G4_Si");
    auto* scintMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

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

    // Thin lithium target. Beam axis = Z.
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

    DetectorLogicals det;
    det.world = logicWorld;

    auto* solidScint = new G4Box("Scintillator",
                                 det.scintHalfThickness,
                                 det.scintHalfWidth,
                                 det.scintHalfHeight);

    det.scintillator = new G4LogicalVolume(solidScint,
                                           scintMat,
                                           "ScintillatorLV");

    auto* solidSilicon = new G4Box("SiliconStrip",
                                   det.siHalfThickness,
                                   det.siHalfWidth,
                                   det.siHalfHeight);

    det.silicon = new G4LogicalVolume(solidSilicon,
                                      siliconMat,
                                      "SiliconStripLV");

    auto* scintVis = new G4VisAttributes(G4Colour(0.1, 0.4, 1.0, 0.28));
    scintVis->SetForceSolid(true);
    det.scintillator->SetVisAttributes(scintVis);

    auto* siliconVis = new G4VisAttributes(G4Colour(0.55, 0.15, 1.0, 0.75));
    siliconVis->SetForceSolid(true);
    det.silicon->SetVisAttributes(siliconVis);

    G4cout << "[DetectorConstruction] Constructing geometry: " << fGeometryName << G4endl;

    switch (fGeometry)
    {
        case GeometryType::Current:
            BuildCurrentAtomkiLikeGeometry(det);
            break;
        case GeometryType::Coverage2Pi:
            BuildTwoPiCoverageGeometry(det);
            break;
        case GeometryType::Coverage4Pi:
            BuildFourPiCoverageGeometry(det);
            break;
        case GeometryType::PadPlane:
            BuildForwardPadPlaneGeometry(det);
            break;
    }

    return physWorld;
}
