#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const G4String& inputFile,
                                           const G4String& outputFile)
    : fInputFile(inputFile),
      fOutputFile(outputFile)
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

    SetUserAction(new PrimaryGeneratorAction(eventAction, fInputFile));
    SetUserAction(new SteppingAction(eventAction));
}
