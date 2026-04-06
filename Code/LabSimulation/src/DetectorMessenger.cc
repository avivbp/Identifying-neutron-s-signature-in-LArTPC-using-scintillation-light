#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det, Run* run)
    : G4UImessenger(), fDetector(det), fRun(run) {


    fDetDir = new G4UIdirectory("/det/");
  fDetDir->SetGuidance("Detector configuration");

  fTopPMTCmd = new G4UIcmdWithAnInteger("/det/setTopPMTs", this);
  fTopPMTCmd->SetGuidance("Active PMTs on top (0..4)");
  fTopPMTCmd->SetParameterName("nTopPMT",false);
  fTopPMTCmd->SetRange("nTopPMT>=0 && nTopPMT<=4");
  fTopPMTCmd->AvailableForStates(G4State_Idle);

  fBotPMTCmd = new G4UIcmdWithAnInteger("/det/setBotPMTs", this);
  fBotPMTCmd->SetGuidance("Active PMTs on bottom (0..4)");
  fBotPMTCmd->SetParameterName("nBotPMT",false);
  fBotPMTCmd->SetRange("nBotPMT>=0 && nBotPMT<=4");
  fBotPMTCmd->AvailableForStates(G4State_Idle);

  fSiPMRowsCmd = new G4UIcmdWithAnInteger("/det/setSiPMRows", this);
  fSiPMRowsCmd->SetGuidance("Active SiPM rows (0..16)");
  fSiPMRowsCmd->SetParameterName("nSiPM",false);
  fSiPMRowsCmd->SetRange("nSiPM>=0 && nSiPM<=16");
  fSiPMRowsCmd->AvailableForStates(G4State_Idle);

  fPDETopPMTCmd = new G4UIcmdWithADouble("/det/setPDETopPMT", this);
  fPDETopPMTCmd->SetGuidance("Top PMT PDE (0..1)");
  fPDETopPMTCmd->SetParameterName("TopPMTPDE",false);
  fPDETopPMTCmd->SetRange("TopPMTPDE>=0.0 && TopPMTPDE<=1.0");
  fPDETopPMTCmd->AvailableForStates(G4State_Idle);

  fPDEBotPMTCmd = new G4UIcmdWithADouble("/det/setPDETopPMT", this);
  fPDEBotPMTCmd->SetGuidance("Top PMT PDE (0..1)");
  fPDEBotPMTCmd->SetParameterName("TopPMTPDE",false);
  fPDEBotPMTCmd->SetRange("TopPMTPDE>=0.0 && TopPMTPDE<=1.0");
  fPDEBotPMTCmd->AvailableForStates(G4State_Idle);

  fPDESiPMCmd = new G4UIcmdWithADouble("/det/setPDESiPM", this);
  fPDESiPMCmd->SetGuidance("SiPM PDE (0..1)");
  fPDESiPMCmd->SetParameterName("SiPMPDE",false);
  fPDESiPMCmd->SetRange("SiPMPDE>=0.0 && SiPMPDE<=1.0");
  fPDESiPMCmd->AvailableForStates(G4State_Idle);

  fReflESRCmd = new G4UIcmdWithADouble("/det/setReflESR", this);
  fReflESRCmd->SetGuidance("ESR reflectivity used when a tile is inactive (0..1)");
  fReflESRCmd->SetParameterName("RefleESR", false);
  fReflESRCmd->SetRange("RefleESR>=0.0 && RefleESR<=1.0");
  fReflESRCmd->AvailableForStates(G4State_Idle);

  fReflPMTCmd = new G4UIcmdWithADouble("/det/setReflPMT", this);
  fReflPMTCmd->SetGuidance("PMT patch reflectivity when active (0..1)");
  fReflPMTCmd->SetParameterName("PMTRefle",false);
  fReflPMTCmd->SetRange("PMTRefle>=0.0 && PMTRefle<=1.0");
  fReflPMTCmd->AvailableForStates(G4State_Idle);

  fReflSiPMCmd = new G4UIcmdWithADouble("/det/setReflSiPM", this);
  fReflSiPMCmd->SetGuidance("SiPM patch reflectivity when active (0..1)");
  fReflSiPMCmd->SetParameterName("SiPMRefle",false);
  fReflSiPMCmd->SetRange("SiPMRefle>=0.0 && SiPMRefle<=1.0");
  fReflSiPMCmd->AvailableForStates(G4State_Idle);

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

    fDetSizeCmd = new G4UIcmdWithADoubleAndUnit("/mysim/setDetSize", this);
    fDetSizeCmd->SetGuidance("Set outer LAr cell diameter");
    fDetSizeCmd->SetParameterName("outerDiameter", false);
    fDetSizeCmd->SetRange("outerDiameter > 0.0");
    fDetSizeCmd->SetUnitCategory("Length");
    fDetSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

DetectorMessenger::~DetectorMessenger() {
  delete fBirksCmd;
  delete fAbsLengthCmd;
  delete fTpbEficCmd;
  delete fDetSizeCmd;
  delete fSimDir;
  delete fTopPMTCmd;
  delete fBotPMTCmd;
  delete fSiPMRowsCmd;
  delete fPDETopPMTCmd;
  delete fPDEBotPMTCmd;
  delete fPDESiPMCmd;
  delete fReflESRCmd;
  delete fReflPMTCmd;
  delete fReflSiPMCmd;
  delete fDetDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == fBirksCmd) {
        fDetector->SetBirksConstant(fBirksCmd->GetNewDoubleValue(newValue));
    } else if (command == fAbsLengthCmd) {
        fDetector->SetAbsorptionLength(fAbsLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == fTpbEficCmd) {
        fRun->SetTpbEfficiency(fTpbEficCmd->GetNewDoubleValue(newValue));
    }
      else if (command == fDetSizeCmd)  {
        G4double newLen = fDetSizeCmd->GetNewDoubleValue(newValue);
        //fDetector->SetOuterCellSize(newLen);
        //fDetector->UpdateOuterCellSize(newLen);
        fDetector->outerDiameter = newLen;
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
    if (command == fTopPMTCmd) {
    fDetector->SetTopPMTs(fTopPMTCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fBotPMTCmd) {
    fDetector->SetBotPMTs(fBotPMTCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fSiPMRowsCmd) {
    fDetector->SetSiPMRows(fSiPMRowsCmd->GetNewIntValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fPDETopPMTCmd) {
    fDetector->SetPDETopPMT(fPDETopPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fPDEBotPMTCmd) {
    fDetector->SetPDEBotPMT(fPDEBotPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }

  else if (command == fPDESiPMCmd) {
    fDetector->SetPDESiPM(fPDESiPMCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  //else if (command == fReflESRCmd) {
    //fDetector->SetReflESR(fReflESRCmd->GetNewDoubleValue(newValue));
    //fDetector->ApplySensorConfig();
 // }
  else if (command == fReflPMTCmd) {
    fDetector->SetReflPMT(fReflPMTCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
  else if (command == fReflSiPMCmd) {
    fDetector->SetReflSiPM(fReflSiPMCmd->GetNewDoubleValue(newValue));
    fDetector->ApplySensorConfig();
  }
}
