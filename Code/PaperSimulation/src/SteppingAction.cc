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
#include "G4EventManager.hh"
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
  //G4EventManager* G4RunManager::eventManager = G4RunManager::GetRunManager();
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
  const G4String & volumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  const std::vector< const G4Track *> * secondaries = aStep->GetSecondaryInCurrentStep();
  G4int size  = (int) (secondaries->size());

  G4ThreeVector prePos = aStep->GetPreStepPoint()->GetPosition()/CLHEP::cm;
  G4ThreeVector postPos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
  G4double absTime = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;

  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
            aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  // if firing neutron only care about elastic events, no more than around 2000 scintillation photons are created
  if (fEventAction->numPhotonsInArgon >= 1000){
      //std::cout << " too many photons" << std::endl;
      G4EventManager::GetEventManager()->AbortCurrentEvent();
  }

  // tracking each scintillation photon individually
  G4int checkByHand = 1;
  if (checkByHand){
  if(particleDef == opticalphoton){
      G4Track* track = aStep->GetTrack();
     
      if(track->GetCreatorProcess()->GetProcessName() == "Scintillation"){

          //if (processName != "Transportation" && processName != "OpAbsorption"){
              //std::cout << processName << std::endl;
          //}

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

          if (aStep->GetTrack()->GetCurrentStepNumber() == 1 && volumeName == "Absorber"){
              fEventAction->numPhotonsInArgon += 1;
          }

          // checking how many scintillation photons reach top base of LAr cyllinder
          if (processName == "Transportation" && volumeName == "Absorber" && fabs(postPos.x()) < 3.75 && fabs(postPos.z()) < 3.75 && trackID != fEventAction->lastTracked ){

              //std::cout << "pre step position = " << prePos << std::endl;
              //std::cout << "post step position = " << postPos << std::endl;
              //std::cout << "post step momentum direction = " << aStep->GetPostStepPoint()->GetMomentumDirection() << std::endl;
              // QE = probability of scintillation photons being abosrbed by the PMTS, estimated 15% as written in Creus et al.
              G4double randNum = (double) rand()/RAND_MAX;
              if (randNum < 0.15){
                  // first p.e signal gives t0 for time of flight calculation
                  // Geant doesn't necessarily look at the photons in time order
		  if (fEventAction->tZeroEx > absTime){
                      fEventAction->tZeroEx = absTime;
                      //std::cout << "tZeroEx = " << absTime << std::endl;
                  }

                  // scintillation photon reached top PMT
                  if (postPos.y() >= 2.35){
                      G4double tUp = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                      fEventAction->tUp = tUp;
                      fEventAction->numReachedUp +=1;
                      fEventAction->lastTracked = trackID;
                  }
                  // scintillation photons reached bottom PMT
                  else if (postPos.y() <= -2.35){
                      G4double tDown = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                      fEventAction->tDown = tDown;
                      fEventAction->numReachedDown +=1;
                      fEventAction->lastTracked = trackID;
                  }
                 
                  G4double cTime = fabs(fEventAction->tUp - fEventAction->tDown);
                  if (fEventAction->tUp > 0 && fEventAction->tDown > 0 && cTime < 200){
                      fEventAction->coincidenceTime = cTime;
                      fEventAction->Coincedence = true;
                  }
                  if (fabs(postPos.y()) >= 2.35){
                      //fEventAction->lastTracked = trackID;
                      //std::cout << "optical photon hit location = " << aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm << " cm " << std::endl;
                      //fEventAction->numPE += 1;
                      G4double tAbsorbed = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                      run->csvfi.open("DecayTime.csv",std::ios_base::app);
                      run->csvfi << eventID << "," << tAbsorbed << "," << std::endl;
                      run->csvfi.close();
                  }
              }
              }

              //only go here if a photon was not detected in the liquid scintillator yet
              if (fEventAction->tOneEx > absTime){ 
                  G4ThreeVector liquidScintPos = fDetector->position2/CLHEP::cm;
                  G4double liquidScintY = liquidScintPos.y();
                  G4double liquidScintX = liquidScintPos.x();
                  G4double lsRad = pow(6.35,2);
                  G4double photRad = pow((liquidScintY - postPos.y()),2) + pow((liquidScintPos.z() - postPos.z()),2);
                  //std::cout << "liquid scint position = " << liquidScintPos << std::endl;
                  if (processName == "Transportation" && volumeName == "liquidScintillator"){ 
                      // photon going outside furthest base of liquid scintillator = "into PMT"
                      if (postPos.x() - liquidScintX >= 6.35*cos(CLHEP::pi/2 - fDetector->rotationAngle)){
                          //std::cout << "pos = " << postPos << std::endl;
                          G4double randNum = (double) rand()/RAND_MAX;
                          // assumed QE for liquid scintillator is 25%
                          if (randNum < 0.25){
                              //std::cout << "time of liquid scint signal = " << aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns << " ns" << std::endl;
                              fEventAction->tOneEx = absTime;
                          }
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
  G4int yep = 1;
  if (yep && trackID == 1 || (parentID == 1 && neutron.compare(incomingParticleName) == 0)){

      //if (volumeName == "liquidScintillator2"){
          //G4double prKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double poKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double en_d = prKineticEnergy - poKineticEnergy;
          //fEventAction->enDep += en_d;
          //if (trans.compare(processName) != 0) {
              //const G4Track* secondary = (*secondaries)[0];
              //if (secondary->GetParticleDefinition()->GetParticleName().compare(proton) == 0) {
                   //G4double en_trans = prKineticEnergy - poKineticEnergy;
                   //fEventAction->enDep += en_trans;
              //}
          //}
          //printEventStats(aStep);
      //}

      //if (volumeName == "Absorber"){
          //G4double preKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double ed = preKineticEnergy - postKineticEnergy;
          //fEventAction->eDep += ed;
      //}

      //printEventStats(aStep);
      //G4double stepLength = calcStepLength(aStep);
      //fEventAction->sumStepLength += stepLength;
      //G4ThreeVector pos = aStep->GetTrack()->GetPosition()/CLHEP::cm;
      //fEventAction->csvfile.open("primaryPos.csv",std::ios_base::app);
      //fEventAction->csvfile << eventID << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," << std::endl;
      //fEventAction->csvfile.close();


      if (volumeName.compare(sensitiveLAr) == 0 && trans.compare(processName) != 0){
          fEventAction->tZero = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
      }

      // elastic scatters of primary neutron with sensitive LAr detector
      if (volumeName.compare(sensitiveLAr) == 0 && elastic.compare(processName) == 0){

          fEventAction->numElasticSensitive += 1;
          fEventAction->scatteredElastically = true;
          G4double preKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //const G4Track * secondary = (*secondaries)[0];
          //G4double Beta = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(), secondary->GetMomentumDirection());
          G4double neutronAngle = calcAngle(aStep->GetPreStepPoint()->GetMomentumDirection(),aStep->GetPostStepPoint()->GetMomentumDirection());

          if (fEventAction->numElasticSensitive == 1){
              fEventAction->nucleusRecoilEnergy = preKineticEnergy - postKineticEnergy;
              fEventAction->nScatterAngle = neutronAngle;
          }

          else if (fEventAction->numElasticSensitive > 1){
              fEventAction->nucleusRecoilEnergy += (preKineticEnergy - postKineticEnergy);
          }

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
    
      if (inelastic.compare(processName) == 0 && volumeName.compare(sensitiveLAr) == 0 ){
          fEventAction->scatteredInelastically = 1;
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

