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

G4ThreeVector RotateVectorAroundAxis(const G4ThreeVector& v,
                                     const G4ThreeVector& axis,
                                     G4double angle)
{
    const G4ThreeVector k = axis.unit();

    const G4double c = std::cos(angle);
    const G4double s = std::sin(angle);

    return v * c
         + k.cross(v) * s
         + k * (k.dot(v) * (1.0 - c));
}

void ApplyWorldEulerRotations(G4ThreeVector& xAxis,
                              G4ThreeVector& yAxis,
                              G4ThreeVector& zAxis,
                              G4double rotX,
                              G4double rotY,
                              G4double rotZ)
{
    // World rotations applied after the default XY-plane orientation.
    //
    // rotX: rotation around global/world X axis.
    // rotY: rotation around global/world Y axis.
    // rotZ: rotation around global/world Z axis.
    //
    // Order:
    //
    //   X, then Y, then Z.
    //
    // This is deliberately world-based, not local-based, so that the rotations
    // correspond directly to the geometrical planes:
    //
    //   YZ plane -> rotation around X
    //   XZ plane -> rotation around Y
    //   XY plane -> rotation around Z

    const G4ThreeVector worldX(1.0, 0.0, 0.0);
    const G4ThreeVector worldY(0.0, 1.0, 0.0);
    const G4ThreeVector worldZ(0.0, 0.0, 1.0);

    if (std::abs(rotX) > 0.0)
    {
        xAxis = RotateVectorAroundAxis(xAxis, worldX, rotX).unit();
        yAxis = RotateVectorAroundAxis(yAxis, worldX, rotX).unit();
        zAxis = RotateVectorAroundAxis(zAxis, worldX, rotX).unit();
    }

    if (std::abs(rotY) > 0.0)
    {
        xAxis = RotateVectorAroundAxis(xAxis, worldY, rotY).unit();
        yAxis = RotateVectorAroundAxis(yAxis, worldY, rotY).unit();
        zAxis = RotateVectorAroundAxis(zAxis, worldY, rotY).unit();
    }

    if (std::abs(rotZ) > 0.0)
    {
        xAxis = RotateVectorAroundAxis(xAxis, worldZ, rotZ).unit();
        yAxis = RotateVectorAroundAxis(yAxis, worldZ, rotZ).unit();
        zAxis = RotateVectorAroundAxis(zAxis, worldZ, rotZ).unit();
    }

    // Re-orthonormalize for numerical safety.
    xAxis = xAxis.unit();
    yAxis = (yAxis - xAxis * yAxis.dot(xAxis)).unit();
    zAxis = xAxis.cross(yAxis).unit();
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
    // Generic helper for planar geometries.
    //
    // Detector boxes are defined as:
    //
    //   G4Box(..., halfThickness, halfWidth, halfHeight)
    //
    // Therefore local +X is the detector normal.

    const G4ThreeVector xAxis = direction.unit();

    const G4ThreeVector zRef(0.0, 0.0, 1.0);
    const G4ThreeVector yRef(0.0, 1.0, 0.0);

    const G4ThreeVector ref =
        (std::abs(xAxis.dot(zRef)) < 0.95) ? zRef : yRef;

    const G4ThreeVector yAxis = ref.cross(xAxis).unit();
    const G4ThreeVector zAxis = xAxis.cross(yAxis).unit();

    auto* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);

    return rot;
}

G4RotationMatrix* RotationForSphericalShell(const G4ThreeVector& direction,
                                            G4double rotXDeg = 0.0,
                                            G4double rotYDeg = 0.0,
                                            G4double rotZDeg = 0.0,
                                            G4bool flipInPlane = false)
{
    // Orientation used for 2pi / 4pi.
    //
    // IMPORTANT:
    //
    // Default state before manual rotations:
    //
    //   detector plane = XY plane
    //   detector normal = global +Z
    //
    // Since detector boxes are defined as:
    //
    //   G4Box(..., halfThickness, halfWidth, halfHeight)
    //
    // local +X is the detector normal.
    //
    // Therefore the default axis mapping is:
    //
    //   local +X -> global +Z
    //   local +Y -> global +X
    //   local +Z -> global +Y
    //
    // Then we apply user/world rotations:
    //
    //   rotXDeg around global X
    //   rotYDeg around global Y
    //   rotZDeg around global Z
    //
    // The input direction is kept in the signature because it is useful for
    // consistency with the old code and for possible future logic, but the
    // default orientation no longer uses a radial/spherical basis.

    (void)direction;

    G4ThreeVector xAxis(0.0, 0.0, 1.0);  // local +X = detector normal = global +Z
    G4ThreeVector yAxis(1.0, 0.0, 0.0);  // local +Y = global +X
    G4ThreeVector zAxis(0.0, 1.0, 0.0);  // local +Z = global +Y

    if (flipInPlane)
    {
        yAxis = -yAxis;
        zAxis = -zAxis;
    }

    ApplyWorldEulerRotations(xAxis,
                             yAxis,
                             zAxis,
                             rotXDeg * deg,
                             rotYDeg * deg,
                             rotZDeg * deg);

    auto* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);

    return rot;
}

void PlaceDirectionalTelescope(const DetectorLogicals& det,
                               const G4ThreeVector& direction,
                               G4double siliconRadius,
                               G4int copyNo,
                               const G4String& namePrefix,
                               G4double rotXDeg = 0.0,
                               G4double rotYDeg = 0.0,
                               G4double rotZDeg = 0.0,
                               G4bool flipOrientationOnly = false)
{
    const G4ThreeVector u = direction.unit();

    const G4double scintRadius =
        siliconRadius
        + det.siHalfThickness
        + det.gapSiScint
        + det.scintHalfThickness;

    PlaceTelescopePair(det,
                       siliconRadius * u,
                       scintRadius * u,
                       RotationForSphericalShell(u,
                                                 rotXDeg,
                                                 rotYDeg,
                                                 rotZDeg,
                                                 flipOrientationOnly),
                       copyNo,
                       namePrefix);
}

G4ThreeVector ExtraRotationDegForDetector(G4bool downstreamHemisphere,
                                          G4int ringIndex,
                                          G4int phiIndex,
                                          G4double thetaDeg,
                                          G4double phiDeg,
                                          const G4ThreeVector& detectorPosition,
                                          const G4ThreeVector& siliconPosition,
                                          const G4ThreeVector& scintillatorPosition,
                                          const G4ThreeVector& direction)
{
    // ============================================================
    // ROTACIONES MANUALES SEGÚN LA POSICIÓN REAL DEL DETECTOR
    // ============================================================
    //
    // Default:
    //
    //   cada detector empieza en el plano XY
    //   o sea, perpendicular a Z
    //
    // Después aplicamos:
    //
    //   plano XZ: tan(alpha_xz) = z / x
    //   plano YZ: tan(alpha_yz) = z / y
    //
    // En código:
    //
    //   alpha_xz = atan2(z, x)
    //   alpha_yz = atan2(z, y)
    //
    // Correspondencia:
    //
    //   plano YZ -> rotación alrededor de X
    //   plano XZ -> rotación alrededor de Y
    //
    // La función devuelve:
    //
    //   G4ThreeVector(rotXDeg, rotYDeg, rotZDeg)
    //
    // donde las rotaciones se aplican alrededor de ejes globales/world.
    //
    // ============================================================

    const G4double x = detectorPosition.x();
    const G4double y = detectorPosition.y();
    const G4double z = detectorPosition.z();
    const G4double r =  std::sqrt(std::pow(x,2) + std::pow(y,2));


    const G4double theta = std::atan2(r, z) / deg ;
    const G4double phi   = std::atan2(y, x) / deg ;

    const G4double rotXDeg = 0 ; // 90-phi;
    const G4double rotYDeg = 90-phi;
    const G4double rotZDeg = 270+theta;

    (void)downstreamHemisphere;
    (void)ringIndex;
    (void)phiIndex;
    (void)thetaDeg;
    (void)phiDeg;
    (void)siliconPosition;
    (void)scintillatorPosition;
    (void)direction;

    return G4ThreeVector(rotXDeg, rotYDeg, rotZDeg);
}

G4int BuildHemisphereShell(const DetectorLogicals& det,
                           G4int firstCopyNo,
                           G4double siliconRadius,
                           G4bool downstreamHemisphere,
                           const G4String& namePrefix,
                           G4bool staggerPhi = false,
                           G4bool flipOrientationPattern = false)
{
    struct RingSpec
    {
        G4double thetaFromAxisDeg;
        G4int nPhi;
        G4double phiOffsetDeg;
    };

    const std::vector<RingSpec> rings = {
        // No detector exactly on the beam axis.
        {28.0,  6,  0.0},
        {55.0, 10, 18.0},
        {80.0, 14,  0.0}
    };

    const G4double scintRadius =
        siliconRadius
        + det.siHalfThickness
        + det.gapSiScint
        + det.scintHalfThickness;

    G4int copyNo = firstCopyNo;

    for (G4int iring = 0; iring < static_cast<G4int>(rings.size()); ++iring)
    {
        const auto& ring = rings[iring];

        const G4double thetaDeg =
            downstreamHemisphere
            ? ring.thetaFromAxisDeg
            : 180.0 - ring.thetaFromAxisDeg;

        const G4double theta = thetaDeg * deg;

        const G4double halfStepDeg =
            staggerPhi ? 0.5 * 360.0 / ring.nPhi : 0.0;

        for (G4int i = 0; i < ring.nPhi; ++i)
        {
            const G4double phiDeg =
                ring.phiOffsetDeg
                + halfStepDeg
                + i * 360.0 / ring.nPhi;

            const G4double phi = phiDeg * deg;

            const G4ThreeVector direction(std::sin(theta) * std::cos(phi),
                                          std::sin(theta) * std::sin(phi),
                                          std::cos(theta));

            const G4ThreeVector u = direction.unit();

            const G4ThreeVector siliconPosition =
                siliconRadius * u;

            const G4ThreeVector scintillatorPosition =
                scintRadius * u;

            const G4ThreeVector detectorPosition =
                (siliconPosition + scintillatorPosition) * 0.5;

            const G4ThreeVector extraRotDeg =
                ExtraRotationDegForDetector(downstreamHemisphere,
                                            iring,
                                            i,
                                            thetaDeg,
                                            phiDeg,
                                            detectorPosition,
                                            siliconPosition,
                                            scintillatorPosition,
                                            direction);

            G4bool flipThisDetector = false;

            if (flipOrientationPattern)
            {
                flipThisDetector =
                    ((direction.x() < 0.0) != (direction.y() < 0.0));
            }

            PlaceDirectionalTelescope(det,
                                      direction,
                                      siliconRadius,
                                      copyNo,
                                      namePrefix,
                                      extraRotDeg.x(),
                                      extraRotDeg.y(),
                                      extraRotDeg.z(),
                                      flipThisDetector);

            ++copyNo;
        }
    }

    return copyNo;
}

void BuildCurrentAtomkiLikeGeometry(const DetectorLogicals& det)
{
    // Original simplified ATOMKI-like hexagonal telescope array.
    //
    // Keep this geometry with the old hand-tuned rotations.
    // This reproduces the clean hexagonal layout.

    const G4int nDet = 6;

    const G4double hexSide = 2.0 * det.scintHalfWidth;
    const G4double hexApothem = 0.5 * std::sqrt(3.0) * hexSide;

    const G4double scintRadius =
        1.03 * (hexApothem + det.scintHalfThickness);

    const G4double siliconRadius =
        scintRadius
        - det.scintHalfThickness
        - det.gapSiScint
        - det.siHalfThickness;

    for (G4int i = 0; i < nDet; ++i)
    {
        const G4double phi = i * 360.0 * deg / nDet;

        const G4double xs = 1.1 * scintRadius * std::cos(phi);
        const G4double ys = 1.1 * scintRadius * std::sin(phi);

        const G4double xi = 1.1 * siliconRadius * std::cos(phi);
        const G4double yi = 1.1 * siliconRadius * std::sin(phi);

        G4double rotZ = 0.0;

        if (i == 0) rotZ =   0.0 * deg;
        if (i == 1) rotZ = 120.0 * deg;
        if (i == 2) rotZ =  60.0 * deg;
        if (i == 3) rotZ =   0.0 * deg;
        if (i == 4) rotZ = 120.0 * deg;
        if (i == 5) rotZ =  60.0 * deg;

        auto* rot = new G4RotationMatrix();
        rot->rotateZ(rotZ);

        PlaceTelescopePair(det,
                           G4ThreeVector(xi, yi, 0.0),
                           G4ThreeVector(xs, ys, 0.0),
                           rot,
                           i,
                           "Current");
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

    /*
     * Queremos un plano de pads en XY:
     *
     *   x cambia
     *   y cambia
     *   z = constante
     *
     * Pero tus detectores están definidos como:
     *
     *   G4Box(..., halfThickness, halfWidth, halfHeight)
     *
     * Por tanto:
     *
     *   eje local X = espesor del detector = normal del detector
     *   cara del detector = plano local YZ
     *
     * Para que la cara local YZ quede en el plano global XY,
     * hay que hacer:
     *
     *   local X -> global +Z   si zSilicon > 0
     *   local X -> global -Z   si zSilicon < 0
     *
     * Eso se consigue con una rotación de +/- 90 grados alrededor de Y.
     */

    const G4double zScint =
        zSilicon
        + sign * (det.siHalfThickness
                  + det.gapSiScint
                  + det.scintHalfThickness);

    G4int copyNo = firstCopyNo;

    for (G4int ix = 0; ix < nx; ++ix)
    {
        for (G4int iy = 0; iy < ny; ++iy)
        {
            const G4double x =
                (ix - 0.5 * (nx - 1)) * pitchX;

            const G4double y =
                (iy - 0.5 * (ny - 1)) * pitchY;

            // Central beam hole.
            if (std::abs(x) < 1.e-9 * mm &&
                std::abs(y) < 1.e-9 * mm)
            {
                continue;
            }

            auto* rot = new G4RotationMatrix();

            /*
             * sign = +1:
             *   local X -> global +Z
             *
             * sign = -1:
             *   local X -> global -Z
             */
            rot->rotateY(-sign * 90.0 * deg);

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
    // Semispherical 2pi coverage on the downstream side (+Z).

    BuildHemisphereShell(det,
                         0,
                         22.0 * cm,
                         true,
                         "TwoPiHemisphere",
                         false,
                         false);
}

void BuildFourPiCoverageGeometry(const DetectorLogicals& det)
{
    // 4pi-inspired layout: two hemispheres around the target.
    //
    // Si quieres desplazar la segunda semiesfera en phi, cambia
    // staggerPhi del segundo BuildHemisphereShell a true.

    G4int nextCopyNo = BuildHemisphereShell(det,
                                            0,
                                            22.0 * cm,
                                            true,
                                            "FourPiForwardHemisphere",
                                            false,
                                            false);

    BuildHemisphereShell(det,
                         nextCopyNo,
                         22.0 * cm,
                         false,
                         "FourPiBackwardHemisphere",
                         false,
                         false);
}

void BuildForwardPadPlaneGeometry(const DetectorLogicals& det)
{
    // Plane of pads downstream of the target, perpendicular to beam axis.

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

    if (key == "current" ||
        key == "atomki" ||
        key == "baseline")
    {
        fGeometry = GeometryType::Current;
        fGeometryName = "current";
    }
    else if (key == "2pi" ||
             key == "coverage2pi" ||
             key == "two_pi")
    {
        fGeometry = GeometryType::Coverage2Pi;
        fGeometryName = "2pi";
    }
    else if (key == "4pi" ||
             key == "coverage4pi" ||
             key == "four_pi")
    {
        fGeometry = GeometryType::Coverage4Pi;
        fGeometryName = "4pi";
    }
    else if (key == "padplane" ||
             key == "pad_plane" ||
             key == "pads" ||
             key == "frontpads")
    {
        fGeometry = GeometryType::PadPlane;
        fGeometryName = "padplane";
    }
    else
    {
        G4ExceptionDescription msg;

        msg << "Unknown geometry '" << geometryName
            << "'. Available geometries: "
            << AvailableGeometries();

        G4Exception("DetectorConstruction::SetGeometry",
                    "X17GEOM001",
                    FatalException,
                    msg);
    }

    G4cout << "[DetectorConstruction] Selected geometry: "
           << fGeometryName << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    auto* worldMat =
        nist->FindOrBuildMaterial("G4_AIR");

    auto* targetMat =
        nist->FindOrBuildMaterial("G4_Li");

    auto* siliconMat =
        nist->FindOrBuildMaterial("G4_Si");

    auto* scintMat =
        nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    auto* solidWorld = new G4Box("World",
                                 70.0 * cm,
                                 70.0 * cm,
                                 70.0 * cm);

    auto* logicWorld = new G4LogicalVolume(solidWorld,
                                           worldMat,
                                           "WorldLV");

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

    // Thin lithium target.
    // Beam axis = Z.

    auto* solidTarget = new G4Tubs("LiTarget",
                                   0.0,
                                   5.0 * mm,
                                   0.05 * mm,
                                   0.0,
                                   360.0 * deg);

    auto* logicTarget = new G4LogicalVolume(solidTarget,
                                            targetMat,
                                            "LiTargetLV");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0.0, 0.0, 0.0),
                      logicTarget,
                      "LiTargetPV",
                      logicWorld,
                      false,
                      0,
                      true);

    auto* targetVis =
        new G4VisAttributes(G4Colour(0.9, 0.2, 0.2));

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

    auto* scintVis =
        new G4VisAttributes(G4Colour(0.1, 0.4, 1.0, 0.28));

    scintVis->SetForceSolid(true);
    det.scintillator->SetVisAttributes(scintVis);

    auto* siliconVis =
        new G4VisAttributes(G4Colour(0.55, 0.15, 1.0, 0.75));

    siliconVis->SetForceSolid(true);
    det.silicon->SetVisAttributes(siliconVis);

    G4cout << "[DetectorConstruction] Constructing geometry: "
           << fGeometryName << G4endl;

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