#include "ActionInitialization.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    auto* runAction = new RunAction();
    SetUserAction(runAction);

    auto* eventAction = new EventAction();
    SetUserAction(eventAction);

    SetUserAction(new PrimaryGeneratorAction(eventAction,
                                             "data/data_pair_creation.txt"));

    SetUserAction(new SteppingAction(eventAction));
}