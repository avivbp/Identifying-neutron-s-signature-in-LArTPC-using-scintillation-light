//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm5/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* DET,
                               EventAction* EA)
:G4UserSteppingAction(),fDetector(DET), fEventAction(EA)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int protonNum = 0;
  G4int neutronNum = 0;
  G4int gammaNum = 0;
  G4int Z = 0;
  G4int A = 0;
  std::string proton("proton");
  std::string neutron("neutron");
  std::string gamma("gamma");
  std::string trans("Transportation");
  std::string elastic("hadElastic");
  std::string inelastic("neutronInelastic");
  G4bool elasticc = false;
  const G4String & processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  
  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION
  if( trackID == 1 && (trans.compare(processName) != 0)){

      // CHECKING IF PRIMARY INTERACTED SEVERAL TIMES
      fEventAction->numInteractionsPP();
      if (fEventAction->numInteractions > 1){
         std::cout << "number of interactions of primary = " << fEventAction->numInteractions << std::endl;
      }

      run->addNumInteractions();
     // std::cout << "process name: " << processName << std::endl;
      if(elastic.compare(processName) == 0){
          run->addNumElastic();
          elasticc = true;
      }

      else if (inelastic.compare(processName) == 0){
          run->addNumInelastic();
          const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
          G4int size  = (int) (secondaries->size());
          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              
              if (proton.compare(particleName) == 0){
                  protonNum+=1;
              }
              else if (neutron.compare(particleName) == 0){
                  neutronNum+=1;
              }
              else if (gamma.compare(particleName) == 0){
                  gammaNum+=1;
              }

              // CHECKING FOR THE HEAVY NUCLEAS AMONG SECONDARIES
              G4int charge = secondary->GetParticleDefinition()->GetAtomicNumber();
              if (charge > 4){
                  Z = charge;
                  A = secondary->GetParticleDefinition()->GetAtomicMass();
              }
          }
          
      }

      // IF HAVEN'T WROTE TO FILE YET, WRITE DATA TO FILE, ELSE IGNORE MORE THAN ONE INTERACTION BY PRIMARY
      if (!(fEventAction->wroteToFile)){
      run-> csvf.open("reactions.csv",std::ios_base::app);
      // run-> csvf << neutronNum << "," << protonNum << "," << gammaNum << "," << Z << "," << A << "," << std::endl;
      if (elasticc){
          run->csvf << "elastic" << "," << std::endl;
      }
      else if(neutronNum == 1 && protonNum == 0){
          run->csvf << "nnprime" << "," << std::endl;
      }
      else if(neutronNum == 1 && protonNum == 1){
          run->csvf << "nnprimeP" << "," << std::endl;
      }
      else if(neutronNum == 2 && protonNum == 0){
          run->csvf << "nTwoN" << "," << std::endl;
      }
      else if(neutronNum == 3 && protonNum == 0){
          run->csvf << "nThreeN" << "," << std::endl;
      }
      else if(neutronNum == 0 && protonNum == 1){
          run->csvf << "np" << "," << std::endl;
      }
      else if(neutronNum == 0 && protonNum == 0 && gammaNum > 0){
         // std::cout << "A of heavy nucleus : " << A << "Z of heavy nucleus : " << Z << "number of gammas: " << gammaNum << std::endl;
          if (A == 37 && Z == 16){
              run->csvf << "alphaKO" << "," << std::endl;
          }
          else{
              run->csvf << "nNGamma" << "," << std::endl;
          }
      }
      else if(neutronNum == 4 && protonNum == 0){
          run->csvf << "nFourN" << "," << std::endl;
      }
      else if(neutronNum == 3 && protonNum == 1){
          run->csvf << "nThreeNP" << "," << std::endl;
      }
      else if(neutronNum == 2 && protonNum == 1){
          run->csvf << "nTwoNP" << "," << std::endl;
      }
      else{
          std::cout << "A of heavy nucleus : " << A << "Z of heavy nucleus : " << Z << "number of gammas: " << gammaNum << "number of neutrons : " << neutronNum << "number of protons: " << protonNum <<  std::endl;
          run->csvf << "other" << "," << std::endl;
      }
      run-> csvf.close();
      fEventAction->setWroteToFile(true);
      }

  }

  if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() 
    != fDetector->GetAbsorber()) return;
    
  fEventAction->AddEnergy (aStep->GetTotalEnergyDeposit());
   
  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  if (charge != 0.) { 
    fEventAction->AddTrakLenCharg(aStep->GetStepLength());
    fEventAction->CountStepsCharg();
  } else {
    fEventAction->AddTrakLenNeutr(aStep->GetStepLength());
    fEventAction->CountStepsNeutr();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

