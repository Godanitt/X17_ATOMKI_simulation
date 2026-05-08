#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const G4String& inputFile,
                                           const G4String& outputFile,
                                           const G4String& generationMode)
    : fInputFile(inputFile),
      fOutputFile(outputFile),
      fGenerationMode(generationMode)
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction(fOutputFile));
}

void ActionInitialization::Build() const
{
    SetUserAction(new RunAction(fOutputFile));

    auto* eventAction = new EventAction();
    SetUserAction(eventAction);

    SetUserAction(new PrimaryGeneratorAction(eventAction, fInputFile, fGenerationMode));
    SetUserAction(new SteppingAction(eventAction));
}
