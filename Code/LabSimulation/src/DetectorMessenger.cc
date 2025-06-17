#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det, Run* run)
    : fDetector(det), fRun(run) {
    fSimDir = new G4UIdirectory("/mysim/");
    fSimDir->SetGuidance("Commands to control simulation parameters");

    fBirksCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setBirks", this);
    fBirksCmd->SetGuidance("Set Birks constant");
    fBirksCmd->SetParameterName("birks", false);
    fBirksCmd->SetRange("birks > 0.0");
    fBirksCmd->SetUnitCategory("Length");
    fBirksCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fAbsLengthCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setAbsLength", this);
    fAbsLengthCmd->SetGuidance("Set absorption length");
    fAbsLengthCmd->SetParameterName("absLength", false);
    fAbsLengthCmd->SetRange("absLength > 0.0");
    fAbsLengthCmd->SetUnitCategory("Length");
    fAbsLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fTpbEficCmd = new G4UIcmdWithADouble("/mysim/setTpbEfic", this);
    fTpbEficCmd->SetGuidance("Set TPB efficiency");
    fTpbEficCmd->SetParameterName("tpbEfic", false);
    fTpbEficCmd->SetRange("tpbEfic >= 0.0 && tpbEfic <= 1.0");
    fTpbEficCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
    delete fBirksCmd;
    delete fAbsLengthCmd;
    delete fTpbEficCmd;
    delete fSimDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fBirksCmd) {
        fDetector->SetBirksConstant(fBirksCmd->GetNewDoubleValue(newValue));
    } else if (command == fAbsLengthCmd) {
        fDetector->SetAbsorptionLength(fAbsLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == fTpbEficCmd) {
        fRun->SetTpbEfficiency(fTpbEficCmd->GetNewDoubleValue(newValue));
    }
}

