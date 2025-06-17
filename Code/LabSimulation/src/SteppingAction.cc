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
#include "G4EventManager.hh"
#include "G4AtomicShells.hh"
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
  G4int parentID = aStep->GetTrack()->GetParentID();
  G4double gammaEDep = 0;
  G4bool nPrimeThisStep = false;
  G4double emittedNeutronAngle = 0.0;
  std::string positron("e+");
  std::string proton("proton");
  std::string neutron("neutron");
  std::string gamma("gamma");
  std::string phot("phot");
  std::string trans("Transportation");
  std::string elastic("hadElastic");
  std::string inelastic("neutronInelastic");
  std::string capture("nCapture");
  std::string scintillator("liquidScintillator");
  std::string sensitiveLAr("Absorber");
  std::string Absorber("biggerAbsorber");

  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  const G4String & processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  const G4String & incomingParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  const G4String & volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
  G4int size  = (int) (secondaries->size());

  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
            aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  G4double pre_Ekin = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
  G4double post_Ekin = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
  G4double postStepTime = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
  G4double preStepTime = aStep->GetPreStepPoint()->GetGlobalTime()/CLHEP::ns;

  const G4ThreeVector & initMomentum = aStep->GetPreStepPoint()->GetMomentumDirection();
  const G4ThreeVector & finMomentum = aStep->GetPostStepPoint()->GetMomentumDirection();

  // tracking each scintillation photon individually
  G4int checkByHand = 0;
  if (checkByHand){
  if(particleDef == opticalphoton){
      G4Track* track = aStep->GetTrack();
     
      if(track->GetCreatorProcess()->GetProcessName() == "Scintillation" ){

          if (aStep->GetTrack()->GetCurrentStepNumber() == 1 && volumeName == "Absorber"){
              fEventAction->numPhotons+=1;

             //if (fEventAction->numPhotons >= 20000){
                  //G4EventManager::GetEventManager()->AbortCurrentEvent();
             //} 
          }   

          // checking how many scintillation photons reach PMTs
          if (processName == "Transportation" && volumeName == "Absorber" && trackID != fEventAction->lastTracked){

              G4double randNum = (double) rand()/RAND_MAX;
              // X-ARAPUCA total efficiency ~ 2% : https://arxiv.org/pdf/2405.12014
              // to save time lowered scint yield times 0.02 instead of QE, change later
              G4double QE = 1.02;

              G4bool reachedPDS = withinArapuca(aStep);
              if (reachedPDS && randNum < QE){

                  fEventAction->tUp = postStepTime;
                  fEventAction->numPE += 1;
                  fEventAction->lastTracked = trackID;
                  if (fEventAction->tZeroEx > preStepTime){
                      fEventAction->tZeroEx = preStepTime;
                  }
              }

          }

      }
  }
  }

  G4int hmm = 1;
  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (hmm && incomingParticleName == "neutron"){

      if (volumeName == "innerLayer" && trans.compare(processName) != 0){
          fEventAction->innerLayerScatter += 1;
          fEventAction->innerLayerEDep += aStep->GetTotalEnergyDeposit()/CLHEP::keV;
      }

      if (volumeName == "outerCell" && trans.compare(processName) != 0){
          fEventAction->ExtScatter += 1;
          fEventAction->outerCellScatterAngle = calcAngle(initMomentum,finMomentum);
          fEventAction->outerCellEDep += aStep->GetTotalEnergyDeposit()/CLHEP::keV;
      }

      if ((volumeName == "outerLayerOne" || volumeName == "outerLayerTwo") && trans.compare(processName) != 0){
          fEventAction->CryoScatter += 1;
          fEventAction->CryoScatterAngle = calcAngle(initMomentum,finMomentum);
          fEventAction->CryoEDep += aStep->GetTotalEnergyDeposit()/CLHEP::keV;
      }


      if (volumeName == "innerCell" && trans.compare(processName) != 0 && fEventAction->tZero > postStepTime){
          fEventAction->tZero = postStepTime;
      }

      if (volumeName == "innerCell" && elastic.compare(processName) == 0){
          fEventAction->numElasticSensitive += 1;
          fEventAction->eDep += (pre_Ekin - post_Ekin);
      }

      if (volumeName == "innerCell" && inelastic.compare(processName) == 0){
          fEventAction->numInelasticSensitive += 1;
      }

      if (inelastic.compare(processName) == 0){ 
          fEventAction->numInelastic += 1;
          //fEventAction->aborted = true;
          //G4EventManager::GetEventManager()->AbortCurrentEvent();
          //run->addNumInelastic();
 
          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              if (neutron.compare(particleName) == 0){
                  fEventAction->secondaryNeutron = true;
                  //emittedNeutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),secondary->GetMomentumDirection());
                  // IF THE PRIMARY NEUTRON PEFFORMED (n,n') LABEL BOOL AS true
                  //fEventAction->setNPrime(true);

                  //G4double energy = secondary->GetKineticEnergy()/CLHEP::MeV; 
                  //run->csvfil.open("nPrimeEnergies.csv",std::ios_base::app);
                  //run->csvfil << eventID << "," << energy << "," << std::endl;
                  //run->csvfil.close();
                 // std::cout << "n' energy = " << energy << " MeV" << std::endl;
             } 
          }

          // only if a neutron came out there was (n,n'), could also be inelastic (n,p) or (n,alpha)
          if (fEventAction->secondaryNeutron){
              nPrimeThisStep = true;
              fEventAction->scatteredInelastically = 1;
          }
       }
       if (processName == "nCapture"){
           fEventAction->nCapture = true;
       }

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

