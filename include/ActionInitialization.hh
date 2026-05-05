#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4String.hh"
#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization
{
public:
    ActionInitialization(const G4String& inputFile,
                         const G4String& outputFile);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4String fInputFile;
    G4String fOutputFile;
};

#endif
