#ifndef X17_EVENT_INFORMATION_HH
#define X17_EVENT_INFORMATION_HH 1

#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"

class X17EventInformation : public G4VUserEventInformation
{
public:
    X17EventInformation(double thetaEE_deg,
                        double thetaEm_deg,
                        double kinEm_MeV,
                        double thetaEp_deg,
                        double kinEp_MeV,
                        const G4ThreeVector& dirEm,
                        const G4ThreeVector& dirEp)
        : fThetaEE(thetaEE_deg),
          fThetaEm(thetaEm_deg),
          fKinEm(kinEm_MeV),
          fThetaEp(thetaEp_deg),
          fKinEp(kinEp_MeV),
          fDirEm(dirEm.unit()),
          fDirEp(dirEp.unit())
    {}

    ~X17EventInformation() override = default;

    void Print() const override {}

    double thetaEE() const { return fThetaEE; }
    double thetaEm() const { return fThetaEm; }
    double thetaEp() const { return fThetaEp; }
    double kinEm() const { return fKinEm; }
    double kinEp() const { return fKinEp; }
    const G4ThreeVector& dirEm() const { return fDirEm; }
    const G4ThreeVector& dirEp() const { return fDirEp; }

private:
    double fThetaEE;
    double fThetaEm;
    double fKinEm;
    double fThetaEp;
    double fKinEp;
    G4ThreeVector fDirEm;
    G4ThreeVector fDirEp;
};

#endif
