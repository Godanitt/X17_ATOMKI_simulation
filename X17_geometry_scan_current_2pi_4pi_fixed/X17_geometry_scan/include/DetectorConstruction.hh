#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"

class G4GenericMessenger;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    explicit DetectorConstruction(const G4String& geometryName = "current");
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    void SetGeometry(const G4String& geometryName);
    const G4String& GetGeometryName() const { return fGeometryName; }

    static G4String AvailableGeometries();

private:
    enum class GeometryType {
        Current,
        Coverage2Pi,
        Coverage4Pi,
        PadPlane
    };

    void DefineCommands();

    GeometryType fGeometry = GeometryType::Current;
    G4String fGeometryName = "current";
    G4GenericMessenger* fMessenger = nullptr;
};

#endif
