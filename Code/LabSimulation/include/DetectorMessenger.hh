#ifndef DetectorMessenger_HH
#define DetectorMessenger_HH

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

class DetectorConstruction;
class Run;

class DetectorMessenger : public G4UImessenger {
public:
    DetectorMessenger(DetectorConstruction* det, Run* run);
    ~DetectorMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    DetectorConstruction* fDetector;
    Run* fRun;

    G4UIdirectory* fSimDir;
    G4UIcmdWithADoubleAndUnit* fBirksCmd;
    G4UIcmdWithADoubleAndUnit* fAbsLengthCmd;
    G4UIcmdWithADouble* fTpbEficCmd;
};

#endif

