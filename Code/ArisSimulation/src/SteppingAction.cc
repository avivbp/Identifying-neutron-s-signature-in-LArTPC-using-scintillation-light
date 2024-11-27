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

  //printEventStats(aStep);

  //std::cout << "in Stepping" << std::endl;
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
  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
  G4double absTime = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;

  // tracking each scintillation photon individually
  G4int checkByHand = 1;
  if (checkByHand){
  if(particleDef == opticalphoton){
      //std::cout << "hi im an optical photon" << std::endl;
      G4Track* track = aStep->GetTrack();
     
      // if firing neutron only care about elastic events, no more than around 2000 scintillation photons are created
      if(track->GetCreatorProcess()->GetProcessName() == "Scintillation" ){

          if (aStep->GetTrack()->GetCurrentStepNumber() == 1){
              fEventAction->numPhotons+=1;
          }   

          // checking how many scintillation photons reach PMTs
          G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
          if (processName == "Transportation" && volumeName.compare(sensitiveLAr) == 0 && trackID != fEventAction->lastTracked){

              if (fabs(pos.y()) != 3.8){
                  std::cout << "evt number " << eventID << std::endl;
                  std::cout << "position isnt on bases for some reason, " << pos << std::endl;
              }

              // QE = probability of scintillation photons being abosrbed by the PMTS
              G4double randNum = (double) rand()/RAND_MAX;
              G4double topQE = 0.266;
              G4double bottomQE = 0.339;
              // efficiency of the TPB converting from 128 nm to 420 nm, see table 1 in https://dspace.mit.edu/bitstream/handle/1721.1/121376/1711.01230.pdf?sequence=1&isAllowed=y
              G4double tpbEfic = 0.44*1.22;
              topQE = topQE * tpbEfic;
              bottomQE = bottomQE * tpbEfic;

              G4double bottomPMTRad = 3.81*CLHEP::cm;

              //effective area of 7 2.05cmX2.05cm PMTs
              G4double topPMTLength = 5.8835*CLHEP::cm;
              G4double topPMTHeight = 5.*CLHEP::cm;

              // scintillation photon reached top PMT array, which have 26.6% QE
              if (pos.y() >= 3.8 && fabs(pos.x()) < topPMTLength/2 && fabs(pos.z()) < topPMTHeight/2 && randNum < topQE){
                  G4double tUp = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                  //run->csvfi.open("absTime.csv",std::ios_base::app);
                  //run->csvfi << eventID << "," << tUp << "," << std::endl;
                  //run->csvfi.close();

                  fEventAction->tUp = tUp;
                  fEventAction->numReachedUp +=1;
                  fEventAction->lastTracked = trackID;
                  if (fEventAction->tZeroEx > absTime){
                      //std::cout << absTime << " ns" << std::endl;
                      fEventAction->tZeroEx = absTime;
                  }
              }

              // scintillation photons reached bottom PMT, which has 33.9% QE
              else if (pos.y() <= -3.8 && (pow(pos.x(),2) + pow(pos.z(),2)) <= pow(bottomPMTRad,2) && randNum < bottomQE){
                  G4double tDown = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                  fEventAction->tDown = tDown;
                  fEventAction->numReachedDown +=1;
                  fEventAction->lastTracked = trackID;
                  if (fEventAction->tZeroEx > absTime){
                      //std::cout << absTime << " ns" << std::endl;                     
                      fEventAction->tZeroEx = absTime;
                  }
              }

              // checking for coincidence between PMTs to make sure it is a valid event
              //G4double cTime = fabs(fEventAction->tUp - fEventAction->tDown);
              //if (fEventAction->tUp > 0 && fEventAction->tDown > 0 && cTime < 100){
                  //fEventAction->coincidenceTime = cTime;
                  //fEventAction->Coincedence = true;
              //}
 
              //std::cout << "optical photon hit location = " << aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm << " cm " << std::endl;
              //fEventAction->numPE += 1;
              //G4double tAbsorbed = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
              //run->csvfi.open("photLoc.csv",std::ios_base::app);
              //run->csvfi << eventID << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," << tAbsorbed << "," << std::endl;
              //run->csvfi.close();
          }

          // timing of photon signal in liquid scintillator for tof
          if ((volumeName == "A0" || volumeName == "A1" || volumeName == "A2" || volumeName == "A3" || volumeName == "A4" || volumeName == "A5" || volumeName == "A6" || volumeName == "A7") && fEventAction->tOneEx > absTime){
              //run->csvfi.open("DecayTime.csv",std::ios_base::app);
              //run->csvfi << eventID << "," << absTime << "," << std::endl;
              //run->csvfi.close();

              if (processName == "Transportation"){
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

  G4int hmm = 1;
  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (hmm && (trackID == 1)){ 

      if (volumeName == "Absorber" && trans.compare(processName) != 0){
          fEventAction->tZero = absTime;
      }

      if (volumeName.compare(sensitiveLAr) == 0 && elastic.compare(processName) == 0){
          fEventAction->numElasticSensitive += 1;
      }


      if (volumeName == "A0" || volumeName == "A1" || volumeName == "A2" || volumeName == "A3" || volumeName == "A4" || volumeName == "A5" || volumeName == "A6" || volumeName == "A7"){

          //std::cout << "eventID = " << eventID << ", neutron hit " << volumeName << " ND" << std::endl;
          //printEventStats(aStep);
          G4double prKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double poKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double en_d = prKineticEnergy - poKineticEnergy;
          //fEventAction->enDep += en_d;

          if (trans.compare(processName) != 0) {
              const G4Track* secondary = (*secondaries)[0];
              if (secondary->GetParticleDefinition()->GetParticleName().compare(proton) == 0) {
                   G4double en_trans = prKineticEnergy - poKineticEnergy;
                   fEventAction->enDep += en_trans;
                   fEventAction->detector = volumeName;
                   fEventAction->tOne = absTime;
                   //std::cout << volumeName << std::endl;
              }
          }
          //printEventStats(aStep);
      }


      // being an experimentalist, not using truth values

      //if (processName.compare(trans) != 0 && volumeName.compare(sensitiveLAr) == 0){
          //fEventAction->tZeroEx = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
          //const G4ThreeVector & pos = aStep->GetPostStepPoint()->GetPosition();
          //fEventAction->pos0 = pos;
      //}

      //if (volumeName.compare(Absorber) == 0 && trans.compare(processName) != 0){
          //fEventAction->scatteredNotSensitive = 1;
      //}
 
      // start timer when the neutron leaves the sensitive volume of LAr
      //if (processName.compare(trans) == 0 && volumeName.compare(sensitiveLAr) == 0) {
          //fEventAction->tZero = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
      //}
    
      if (inelastic.compare(processName) == 0 && (volumeName.compare(sensitiveLAr) == 0 || volumeName.compare(Absorber) == 0)){
          //run->addNumInelastic();
 
          G4bool secondaryNeutron = false;       
          for(int i = 0; i < size; i++){
              const G4Track* secondary = (*secondaries)[i];
              const G4String & particleName = secondary->GetParticleDefinition()->GetParticleName();
              if (neutron.compare(particleName) == 0){
                  secondaryNeutron = true;
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

      // if a neutron has an interaction with the liquid scintillator
      if (!fEventAction->detected && trans.compare(processName) != 0 && (volumeName.compare(scintillator) == 0)) {
          //std::cout << "interaction in liquid scintillator" << std::endl;
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

  //std::cout << "in end of stepping" << std::endl;

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

