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

 // printEventStats(aStep);

  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int parentID = aStep->GetTrack()->GetParentID();
  G4int totalNeutrons = 0;
  G4int totalProtons = 0;
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
  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  const G4String & processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  const G4String & incomingParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
 

  // CHECKING THE NEUTRONS' KINETIC ENERGY AS THEY EXIT THE DETECTOR
  if (neutron.compare(incomingParticleName) == 0 && trans.compare(processName) == 0 && aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber()){
     // std::cout << aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV << std::endl;
     run->csvfi.open("neutronEscapeEnergy.csv",std::ios_base::app);
     run->csvfi << eventID << "," << aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV << "," << fEventAction->performedNPrime << "," << std::endl;

     if (fEventAction->performedNPrime){
         run->csvfi << eventID << "," << 0 << "," << 0 << "," << std::endl;
     }

     run->csvfi.close();
     fEventAction->setWroteToFile(true);
  }


  //CHECKING FOR ALL NEUTRON,PROTONS AND PHOTON CREATED IN AN EVENT
  if (trans.compare(processName) != 0){
     // if ((eventID % 1000) == 0){
       //   printEventStats(aStep);
      //}
 
      const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
      G4int size  = (int) (secondaries->size());

      if (gamma.compare(incomingParticleName) == 0){

         // gammaEDep = aStep->GetTotalEnergyDeposit()/CLHEP::keV;
          gammaEDep = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV - aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;

          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();

              if (gamma.compare(particleName) == 0){
                   G4double kineticEnergy = secondary->GetKineticEnergy()/CLHEP::keV;
                   gammaEDep -= kineticEnergy;
              }
          }
          //printEventStats(aStep);
         // std::cout << "gammaEDep = " << gammaEDep << std::endl;

         // std::cout << "gammaEDep = " << gammaEDep << "," << aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV - aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV << std::endl;
          if (aStep->GetTrack()->GetParentID() == 1){
              fEventAction->addDeposit(gammaEDep);

              run->csvfi.open("photonEDIst.csv",std::ios_base::app);
              run->csvfi << eventID << "," << gammaEDep << "," << std::endl;
              run->csvfi.close();
          }
          fEventAction->addEDep(gammaEDep);
       }
          
      for(int i = 0; i < size; i++){
           const G4Track* secondary = (*secondaries)[i];
           const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
           
           if (proton.compare(particleName) == 0 && secondary->GetKineticEnergy()/CLHEP::MeV > 0.5){
               totalProtons += 1;
           }
           else if (neutron.compare(particleName) == 0){
               totalNeutrons += 1;
           }
           else if (gamma.compare(particleName) == 0){
               G4double kineticEnergy = secondary->GetKineticEnergy()/CLHEP::MeV;
               G4double time = aStep->GetPostStepPoint()->GetGlobalTime();

               run->csvfi.open("gammaEDist.csv",std::ios_base::app);
               run->csvfi << eventID << "," << kineticEnergy << "," << time << "," << std::endl;
               run->csvfi.close();

           }
      }



      //CHECKING IF THE PARTICLE THAT CREATED THE SECONDARIES IS A PROTON OR NEUTRON AND SUMMING THE AMOUNT OF CREATED PARTICLES ACCORDINGLY
      if (proton.compare(incomingParticleName) == 0 && totalProtons != 0){
          fEventAction->addProtonNum(totalProtons - 1);
          fEventAction->addNeutronNum(totalNeutrons);
      }

      else if (neutron.compare(incomingParticleName) == 0 && totalNeutrons != 0){
          fEventAction->addProtonNum(totalProtons);
          fEventAction->addNeutronNum(totalNeutrons - 1);
      }

      else {
          fEventAction->addProtonNum(totalProtons);
          fEventAction->addNeutronNum(totalNeutrons);
      }

      //if (neutron.compare(incomingParticleName) == 0 && totalNeutrons == 0 && elastic.compare(processName) != 0 ){
          //printEventStats(initialParticleName,preEnergy,postEnergy,aStep);
      //}
  }



  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (trackID == 1 || (parentID == 1 && neutron.compare(incomingParticleName) == 0) ){
 
      if (aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV == 2.5 && elastic.compare(processName) == 0){

          G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV;
          const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
          const G4Track * secondary = (*secondaries)[0];
          G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
          G4double neutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());

          run->csvfi.open("elasticEnergies.csv",std::ios_base::app);
          run->csvfi << eventID << "," << postKineticEnergy << "," << Beta << "," << neutronAngle << "," << std::endl;
          run->csvfi.close();
      }
     
      if (elastic.compare(processName) == 0){
          run->addNumElastic();
          fEventAction->elasticPP();
      }

      if (inelastic.compare(processName) == 0){
          run->addNumInelastic();
          const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
          G4int size  = (int) (secondaries->size());
 
          G4bool secondaryNeutron = false;       
          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              if (neutron.compare(particleName) == 0){
                  secondaryNeutron = true;
                  emittedNeutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),secondary->GetMomentumDirection());
                  // IF THE PRIMARY NEUTRON PEFFORMED (n,n') LABEL BOOL AS true
                  fEventAction->setNPrime(true);

                  G4double energy = secondary->GetKineticEnergy()/CLHEP::MeV; 
                  run->csvfil.open("nPrimeEnergies.csv",std::ios_base::app);
                  run->csvfil << eventID << "," << energy << "," << std::endl;
                  run->csvfil.close();
                 // std::cout << "n' energy = " << energy << " MeV" << std::endl;
              }
          }

          // only if a neutron came out there was (n,n'), could also be inelastic (n,p) or (n,alpha)
          if (secondaryNeutron){
              nPrimeThisStep = true;

              for(int i = 0; i < size; i++){

                  const G4Track* secondary = (*secondaries)[i];
                  if (secondary->GetParticleDefinition()->GetAtomicNumber() > 10){

                      G4double nucleusEnergy = secondary->GetKineticEnergy()/CLHEP::keV;
                      G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
                      G4double neutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());
                      run->csvfi.open("inelasticEnergies.csv",std::ios_base::app);
                      run->csvfi << eventID << "," << nucleusEnergy << "," << Beta << "," << emittedNeutronAngle << "," << std::endl;
                      run->csvfi.close();
                  }
              }
          }
       }

      // ONLY SUM STEPS THAT ARE WITHIN THE lAr ABSORBER
      if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber()){

          // COUNT HOW MANY INTERACTIONS PRIMARY PERFORMED
          if (trans.compare(processName) != 0) {
              fEventAction -> numInteractionsPP();
          }
          
          G4double stepLength = calcStepLength(aStep);
          G4double timeDiff = (aStep->GetPostStepPoint()->GetGlobalTime() - aStep->GetPreStepPoint()->GetGlobalTime())/CLHEP::nanosecond;
          fEventAction->addStepLength(stepLength);

          // IF PRIMARY PERFORMED (n,n') STOP TRACKING ITS PATH LENGTH, START TRACKING SECONDARY NEUTRON PATH LENGTH
              if (nPrimeThisStep){
                  csvfile.open("runStats.csv",std::ios_base::app);
                  csvfile << 0 << "," << 0 << "," << fEventAction->sumStepLengths << "," << fEventAction->numInteractions << "," << fEventAction->numElastic << "," << 1 << "," << std::endl;
                  csvfile.close(); 
                  fEventAction->resetNumInteractions();
                  fEventAction->numElastic = 0;      
                  fEventAction->resetStepLengths();
                  nPrimeThisStep = false;
               }

          run->csvfil.open("interactionStats.csv",std::ios_base::app);
          run->csvfil << stepLength << "," << timeDiff << "," << std::endl;
          run->csvfil.close();
      }

     // printEventStats(aStep);
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

