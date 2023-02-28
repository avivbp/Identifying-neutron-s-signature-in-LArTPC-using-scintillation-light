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

  //printEventStats(aStep);

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
  std::string scintillator("liquidScintillator");
  std::string sensitiveLAr("Absorber");
  std::string Absorber("biggerAbsorber");

  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  const G4String & processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  const G4String & incomingParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
  G4int size  = (int) (secondaries->size());

  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
            aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  //if (particleDef == opticalphoton && processName != "Transportation") {
      //std::cout << aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm << std::endl;
  //}

  //if (size>100 && (*secondaries)[0]->GetDynamicParticle()->GetParticleDefinition() == opticalphoton){
      //std::cout << "number of scint photons in this step = " << size << std::endl;
      //std::cout << "proton pos = " << aStep->GetTrack()->GetPosition()/CLHEP::cm << std::endl;
      //std::cout << "proton energy = " << aStep->GetTrack()->GetKineticEnergy()/CLHEP::MeV << std::endl;
  //}

  //if (aStep->GetTrack()->GetTrackID() > fEventAction->largestTrackID){
      //fEventAction->largestTrackID = aStep->GetTrack()->GetTrackID();
      //fEventAction->numPhotons+=1;
  //}  
 
  // tracking each scintillation photon individually
  G4int checkByHand = 1;
  if (checkByHand){
  if(particleDef == opticalphoton){
      G4Track* track = aStep->GetTrack();
     
      //if(track->GetCreatorProcess()->GetProcessName() != "Scintillation"){
          //std::cout << track->GetCreatorProcess()->GetProcessName() << std::endl;
          //fEventAction->nCher+=1;
      //}

      // if firing neutron only care about elastic events, no more than around 2000 scintillation photons are created
      if(track->GetCreatorProcess()->GetProcessName() == "Scintillation" && fEventAction->numPhotons < 2000){

          if (aStep->GetTrack()->GetCurrentStepNumber() == 1){
              fEventAction->numPhotons+=1;

              //G4double energy = track->GetDynamicParticle()->GetKineticEnergy()/CLHEP::eV;
              //G4double time = track->GetGlobalTime()/CLHEP::ns;
              //G4ThreeVector dir = track->GetDynamicParticle()->GetMomentumDirection();
              //G4ThreeVector h2 = G4ThreeVector(dir.x(),0.,dir.z()) / sqrt(dir.x()*dir.x()+dir.z()*dir.z());
              //G4double theta = calcAngle(G4ThreeVector(0.,1.,0.),dir);
              //G4double phi = calcAngle(G4ThreeVector(1.,0.,0.),h2);
              //G4ThreeVector pos = track->GetPosition()/CLHEP::cm;
              //G4ThreeVector pos = track->GetVertexPosition()/CLHEP::cm;


              //fEventAction->csvfile.open("scintCheck.csv",std::ios_base::app);
              //fEventAction->csvfile << energy << "," << time << "," << theta << "," << phi << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," << std::endl;
              //fEventAction->csvfile.close();
          }   

          // checking how many scintillation photons reach top base of LAr cyllinder
          G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
          if (processName == "Transportation" && aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Absorber" && fabs(pos.x()) < 3.75 && fabs(pos.z()) < 3.75 && trackID != fEventAction->lastTracked ){

              // QE = probability of scintillation photons being abosrbed by the PMTS, estimated 15% as written in Creus et al.
              G4double randNum = (double) rand()/RAND_MAX;
              if (randNum < 0.15){
                  
                  // scintillation photon reached top PMT
                  if (pos.y() >= 2.35){
                      G4double tUp = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                      fEventAction->tUp = tUp;
                      fEventAction->numReachedUp +=1;
                      fEventAction->lastTracked = trackID;
                  }
                  // scintillation photons reached bottom PMT
                  else if (pos.y() <= -2.35){
                      G4double tDown = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                      fEventAction->tDown = tDown;
                      fEventAction->numReachedDown +=1;
                      fEventAction->lastTracked = trackID;
                  }
                 
                  G4double cTime = fabs(fEventAction->tUp - fEventAction->tDown);
                  if (fEventAction->tUp > 0 && fEventAction->tDown > 0 && cTime < 100){
                      fEventAction->coincidenceTime = cTime;
                      fEventAction->Coincedence = true;
                  }
              }
              }
          }
  }
  }

  //if (incomingParticleName == "e-" && aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Absorber"){
      //fEventAction->eDep += (aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV - aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::MeV);
  //}

  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (trackID == 1 || (parentID == 1 && neutron.compare(incomingParticleName) == 0)){

      if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "Absorber"){
          G4double preKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double ed = preKineticEnergy - postKineticEnergy;
          fEventAction->eDep += ed;
      }

      //if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "PhysAl"){
          //run->csvfi.open("gammaE.csv",std::ios_base::app);
          //run->csvfi << eventID << "," << aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV << "," << std::endl;
          //run->csvfi.close();
      //}

      //if (processName == "compt") {
          //G4double erEnergy = (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
          //if (erEnergy > 11.4 && erEnergy < 11.6){ 
              //printEventStats(aStep);
          //}
      //}
  
      //printEventStats(aStep);
      //G4double stepLength = calcStepLength(aStep);
      //fEventAction->sumStepLength += stepLength;
      //G4ThreeVector pos = aStep->GetTrack()->GetPosition()/CLHEP::cm;
      //fEventAction->csvfile.open("primaryPos.csv",std::ios_base::app);
      //fEventAction->csvfile << eventID << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," << std::endl;
      //fEventAction->csvfile.close();

      const G4String & volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

      if (volumeName.compare(sensitiveLAr) == 0 && trans.compare(processName) != 0){
          fEventAction->tZero = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
      }

      // elastic scatters of primary neutron with sensitive LAr detector
      if (volumeName.compare(sensitiveLAr) == 0 && elastic.compare(processName) == 0){

          fEventAction->numElasticSensitive += 1;
          fEventAction->scatteredElastically = true;
          //G4double preKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //const G4Track * secondary = (*secondaries)[0];
          //G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
          //G4double neutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());

          //if (fEventAction->numElasticSensitive == 1){
              //fEventAction->nucleusRecoilEnergy = preKineticEnergy - postKineticEnergy;
              //fEventAction->nScatterAngle = neutronAngle;
          //}

          //else if (fEventAction->numElasticSensitive > 1){
              //fEventAction->nucleusRecoilEnergy += (preKineticEnergy - postKineticEnergy);
          //}

          //run->csvfi.open("elasticEnergies.csv",std::ios_base::app);
          //run->csvfi << eventID << "," << postKineticEnergy << "," << Beta << "," << neutronAngle << "," << std::endl;
          //run->csvfi.close();
      }

      // being an experimentalist, not using truth values

      //if (processName.compare(trans) != 0 && volumeName.compare(sensitiveLAr) == 0){
          //fEventAction->tZeroEx = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
          //const G4ThreeVector & pos = aStep->GetPostStepPoint()->GetPosition();
          //fEventAction->pos0 = pos;
      //}

      if (volumeName.compare(Absorber) == 0 && trans.compare(processName) != 0){
          fEventAction->scatteredNotSensitive = 1;
      }
 
      // start timer when the neutron leaves the sensitive volume of LAr
      //if (processName.compare(trans) == 0 && volumeName.compare(sensitiveLAr) == 0) {
          //fEventAction->tZero = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
      //}
    
      if (inelastic.compare(processName) == 0 && (volumeName.compare(sensitiveLAr) == 0 || volumeName.compare(Absorber) == 0)){
          //run->addNumInelastic();
 
          G4bool secondaryNeutron = false;       
          //for(int i = 0; i < size; i++){
              //const G4Track* secondary = (*secondaries)[i];
              //const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              //if (neutron.compare(particleName) == 0){
                  //secondaryNeutron = true;
                  //emittedNeutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),secondary->GetMomentumDirection());
                  // IF THE PRIMARY NEUTRON PEFFORMED (n,n') LABEL BOOL AS true
                  //fEventAction->setNPrime(true);

                  //G4double energy = secondary->GetKineticEnergy()/CLHEP::MeV; 
                  //run->csvfil.open("nPrimeEnergies.csv",std::ios_base::app);
                  //run->csvfil << eventID << "," << energy << "," << std::endl;
                  //run->csvfil.close();
                 // std::cout << "n' energy = " << energy << " MeV" << std::endl;
              //}
          //}

          // only if a neutron came out there was (n,n'), could also be inelastic (n,p) or (n,alpha)
          if (secondaryNeutron){
              nPrimeThisStep = true;
              fEventAction->scatteredInelastically = 1;

              //for(int i = 0; i < size; i++){

                  //const G4Track* secondary = (*secondaries)[i];
                  //if (secondary->GetParticleDefinition()->GetAtomicNumber() > 10){

                      //G4double nucleusEnergy = secondary->GetKineticEnergy()/CLHEP::keV;
                      //G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
                      //G4double neutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());
                      //run->csvfi.open("inelasticEnergies.csv",std::ios_base::app);
                      //run->csvfi << eventID << "," << nucleusEnergy << "," << Beta << "," << emittedNeutronAngle << "," << std::endl;
                      //run->csvfi.close();

                      //fEventAction->nucleusRecoilEnergy = nucleusEnergy;
                      //fEventAction->nScatterAngle = emittedNeutronAngle;
                  //}
              //}
          }
       }

       if (inelastic.compare(processName) == 0 && volumeName.compare(Absorber) == 0) {
           fEventAction->extInelastic = true;
       }

      // if a neutron has an interaction with the liquid scintillator
      if (!fEventAction->detected && trans.compare(processName) != 0 && (volumeName.compare(scintillator) == 0)) {
          const G4Track* secondary = (*secondaries)[0]; 
          if (secondary->GetParticleDefinition()->GetParticleName().compare(proton) == 0 && secondary->GetKineticEnergy()/CLHEP::keV > 680) {
              //printEventStats(aStep);
              fEventAction->detected = true;
              //fEventAction->pos1 = aStep->GetPostStepPoint()->GetPosition();

              //const G4ThreeVector & initMomentumDirection = G4ThreeVector(1.,0.,0.);
              //const G4ThreeVector & finMomentumDirection = (fEventAction->pos1 - fEventAction->pos0).unit();
              //G4double scatAngle = calcAngle(initMomentumDirection,finMomentumDirection);
              //fEventAction->experimentalRecoilEnergy = calcRecoilEnergy(2450*CLHEP::keV,40,scatAngle);

              // stop timer when neutron deposits energy in the organic scintillator
              fEventAction->tOne = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
          }  
      }	
  }


  if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume() != fDetector->GetAbsorber()) return;
    
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

