#ifndef DetectorMessenger_HH
#define DetectorMessenger_HH

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4RunManager.hh"

class DetectorConstruction;
class Run;

class DetectorMessenger : public G4UImessenger {
public:
    DetectorMessenger(DetectorConstruction* det, Run* run);
    ~DetectorMessenger();

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    DetectorConstruction* fDetector;
    G4UIdirectory* fDetDir = nullptr;
    Run* fRun;

    // integer config
    G4UIcmdWithAnInteger* fTopPMTCmd  = nullptr;
    G4UIcmdWithAnInteger* fBotPMTCmd  = nullptr;
    G4UIcmdWithAnInteger* fSiPMRowsCmd = nullptr;

    // doubles
    G4UIcmdWithADouble* fPDETopPMTCmd    = nullptr;
    G4UIcmdWithADouble* fPDEBotPMTCmd    = nullptr;
    G4UIcmdWithADouble* fPDESiPMCmd   = nullptr;

    G4UIcmdWithADouble* fReflESRCmd   = nullptr;
    G4UIcmdWithADouble* fReflPMTCmd   = nullptr;
    G4UIcmdWithADouble* fReflSiPMCmd  = nullptr;  

    G4UIdirectory* fSimDir;
    G4UIcmdWithADoubleAndUnit* fBirksCmd;
    G4UIcmdWithADoubleAndUnit* fAbsLengthCmd;
    G4UIcmdWithADouble* fTpbEficCmd;
    G4UIcmdWithADoubleAndUnit* fDetSizeCmd;
};

#endif

