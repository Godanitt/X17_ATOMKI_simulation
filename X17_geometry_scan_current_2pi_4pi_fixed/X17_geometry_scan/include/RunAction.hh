#ifndef RunAction_h
#define RunAction_h 1

#include "G4String.hh"
#include "G4UserRunAction.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    explicit RunAction(const G4String& outputFile);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

private:
    G4String fOutputFile;
};

#endif
