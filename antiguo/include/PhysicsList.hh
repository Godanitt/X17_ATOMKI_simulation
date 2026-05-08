#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"

class PhysicsList : public G4VUserPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override = default;

protected:
    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;
};

#endif
