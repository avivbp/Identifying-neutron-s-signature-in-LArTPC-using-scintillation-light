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

  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  const G4ParticleDefinition* particleDef =
            aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition();

  G4double pre_Ekin = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
  G4double post_Ekin = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;

  if(particleDef == opticalphoton && processName.compare(trans) == 0 && volumeName.compare(Absorber) == 0){
      fEventAction->numPE += 1;
  }


  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;

  if (incomingParticleName == "gamma"){
      G4double E_diff = pre_Ekin - post_Ekin;
      fEventAction->eDep =+ E_diff;
  }

  //printEventStats(aStep);
  //G4bool secondary_e = false;
  //if (processName != "Transportation"){
      //for(int i = 0; i < size; i++){ 
          //const G4Track* secondary = (*secondaries)[i];
          //if (secondary->GetParticleDefinition()->GetParticleName() == "e-"){
            //  secondary_e = true;
          //}
      //}
      //if (secondary_e){
          //printEventStats(aStep);
       //   secondary_e = false;
      //}
 // }
  //G4double init_Ek = 59.5;

  //if (incomingParticleName == "neutron" && processName.compare(elastic) == 0){
    //  fEventAction->numElastic += 1;
      //G4double recoilEnergy = (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
      //if (recoilEnergy > 200.1){
        //  printEventStats(aStep);
     // }
     // run->csvfi.open("VitaliyStats.csv",std::ios_base::app);
      //run->csvfi << eventID << "," << recoilEnergy << "," << std::endl;
      //run->csvfi.close();
 // }

  if (pre_Ekin == 59.5 && processName == "phot" && size >= 1){

      //run->csvfi.open("numPhotSec.csv",std::ios_base::app);
      //run->csvfi << eventID << "," << size << "," << std::endl;
      //run->csvfi.close();

      //printEventStats(aStep);
      //std::cout << "number of optical photons created = " << size << std::endl;

      G4int num_shells = G4AtomicShells::GetNumberOfShells(18);
      //std::cout << "number of shells in Ar = " << num_shells << std::endl;
      for (int i = 0; i < 7; i++){
          G4double shellEnergy = G4AtomicShells::GetBindingEnergy(18,i);
          G4int num_electrons = G4AtomicShells::GetNumberOfElectrons(18,i);
          //std::cout << "Ar level" << i << " binding energy = " << shellEnergy/CLHEP::eV << " eV" << ", number of electrons in the shell = " << num_electrons << std::endl;
      }
      //printEventStats(aStep);
      G4double E_diff = pre_Ekin - (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
      //run->csvfi.open("energy_diff.csv",std::ios_base::app);
      //run->csvfi << eventID << "," << E_diff << "," << std::endl;
      //run->csvfi.close();

      if (E_diff < 0.1776){
          //printEventStats(aStep);
          //std::cout << "energy diff between primary gamma and energetic electron = " << E_diff << " keV" << std::endl;
      }
      //for(int i = 0; i < size; i++){
          //const G4Track* secondary = (*secondaries)[i];
         // G4double E_k = secondary->GetKineticEnergy()/CLHEP::keV;
          //std::cout << E_k << std::endl;
        //  init_Ek -= E_k;
      //}
      //run->csvfi.open("energy_diff.csv",std::ios_base::app);
     // run->csvfi << eventID << "," << init_Ek << "," << std::endl;
   //   run->csvfi.close();

      //std::cout << "remaining energy from photoelectric effect = " << init_Ek*1000 << " eV" << std::endl;
      //G4double e_kin = (*secondaries)[0]->GetKineticEnergy()/CLHEP::keV;
      //run->csvfi.open("phot_eKin.csv",std::ios_base::app);
      //run->csvfi << e_kin << "," << std::endl;
      //run->csvfi.close();
  }
  //std::cout << "in stepping3" << std::endl;

  // tracking each scintillation photon individually
  G4int checkByHand = 1;
  if (checkByHand){
  if(particleDef == opticalphoton){
      //std::cout << "hi im an optical photon" << std::endl;
      G4Track* track = aStep->GetTrack();
     
      if(track->GetCreatorProcess()->GetProcessName() == "Scintillation" ){

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

          // checking how many scintillation photons reach PMTs
          G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm;
          if (processName == "Transportation" && volumeName.compare(sensitiveLAr) == 0 && trackID != fEventAction->lastTracked){

              // QE = probability of scintillation photons being abosrbed by the PMTS
              G4double randNum = (double) rand()/RAND_MAX;
              G4double QE = 0.339;
              // efficiency of conversion from 128 nm to 420 nm, see table 1 in https://dspace.mit.edu/bitstream/handle/1721.1/121376/1711.01230.pdf?sequence=1&isAllowed=y
              G4double tpbEfic = 0.4 * 1.22; 
              tpbEfic = 0.44*1.22;
              // total efficiency taking into account the QE of the PMT and the QE of the tpb light conversion
              G4double totEfficiency = QE * tpbEfic;
              // effective radius of 7 3" (diameter) PMTs
              G4double PMTRad = 10.08*CLHEP::cm;

              // scintillation photon reached top PMT array, which have 33.9% QE
              if (pos.y() >= 11.75 && (pow(pos.x(),2) + pow(pos.z(),2)) <= pow(PMTRad,2) && randNum < totEfficiency){

                  G4double tUp = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                  fEventAction->tUp = tUp;
                  fEventAction->numReachedUp +=1;
                  fEventAction->lastTracked = trackID;
              }

              // scintillation photons reached bottom PMT array, which have 33.9% QE
              else if (pos.y() <= -11.75 && (pow(pos.x(),2) + pow(pos.z(),2)) <= pow(PMTRad,2) && randNum < totEfficiency){
                  G4double tDown = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
                  fEventAction->tDown = tDown;
                  fEventAction->numReachedDown +=1;
                  fEventAction->lastTracked = trackID;
              }
              
              // current definition of tUp,tDown isnt right, need to check for the fastest photons, they arent necessarily registered first
              G4double cTime = fabs(fEventAction->tUp - fEventAction->tDown);
              if (fEventAction->tUp > 0 && fEventAction->tDown > 0 && cTime < 100){
                  fEventAction->coincidenceTime = cTime;
                  fEventAction->Coincedence = true;
              }
 
              //std::cout << "optical photon hit location = " << aStep->GetPostStepPoint()->GetPosition()/CLHEP::cm << " cm " << std::endl;
              //fEventAction->numPE += 1;
              //G4double tAbsorbed = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
              //run->csvfi.open("photLoc.csv",std::ios_base::app);
              //run->csvfi << eventID << "," << pos.x() << "," << pos.y() << "," << pos.z() << "," << tAbsorbed << "," << std::endl;
              //run->csvfi.close();
          }
      }
  }
  }

  // checking for single elastic scatter with 45-55 keV recoil nucleus
  if (0){
  if (trackID == 1 && trans.compare(processName) != 0 && processName != "NoProcess" && volumeName == "biggerAbsorber"){
      if (elastic.compare(processName) == 0){
          //printEventStats(aStep);
          const G4Track * secondary = (*secondaries)[0];
          G4double kinEn = secondary->GetKineticEnergy()/CLHEP::keV;
          if ((fEventAction->numElasticSensitive >= 1 || kinEn < 203 || kinEn > 226)){
              fEventAction->fail = 1;
              fEventAction->numElasticSensitive += 1;
              //std::cout << "fail" << std::endl;
              //G4EventManager::GetEventManager()->AbortCurrentEvent();
          }
          else{
              fEventAction->numElasticSensitive += 1;
              fEventAction->tZero = aStep->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
          }
      }
      // any interaction other than elastic - can terminate event
      else{
          if (inelastic.compare(processName) == 0){
              //fEventAction->numInelastic += 1;
          }
          fEventAction->fail = 1;
          //std::cout << "fail" << std::endl;
      }
  }
  }

  G4int hmm = 0;
  G4double psEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
  if (hmm && parentID == 1 && "Ar40" == incomingParticleName && aStep->GetTrack()->GetCreatorProcess()->GetProcessName() == "hadElastic" && psEnergy > 45 && psEnergy < 55){

      if (processName == "ionIoni"){
          //std::cout << "number of photons created = " << size << std::endl;
          run->csvfi.open("elasticPhot.csv",std::ios_base::app);
          run->csvfi << aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV << "," << size << "," << aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV << "," << std::endl;
          run->csvfi.close();
          //std::cout << "post Step nucleus energy = " << aStep->GetPostStepPoint()->GetKineticEnergy() << std::endl;
      }

      //std::cout << aStep->GetTrack()->GetCreatorProcess()->GetProcessName() << std::endl;
  } 

  hmm = 0;
  // PRIMARY PARTICLE UNDERTAKING ANY PROCESS THAT ISNT TRANSPORTATION ON ITS WAY TO THE DETECTOR OR ON ITS WAY OUT OF THE DETECTOR
  if (hmm && (trackID == 1 || (parentID == 1 && neutron.compare(incomingParticleName) == 0))){

      if (volumeName == "liquidScintillator2"){
          G4double prKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double poKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          //G4double en_d = prKineticEnergy - poKineticEnergy;
          //fEventAction->enDep += en_d;
          if (trans.compare(processName) != 0) {
              const G4Track* secondary = (*secondaries)[0];
              if (secondary->GetParticleDefinition()->GetParticleName().compare(proton) == 0) {
                   G4double en_trans = prKineticEnergy - poKineticEnergy;
                   fEventAction->enDep += en_trans;
              }
          }
          //printEventStats(aStep);
      }

      if (volumeName == "Absorber"){
          G4double preKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double postKineticEnergy = aStep->GetPostStepPoint()->GetKineticEnergy()/CLHEP::keV;
          G4double ed = preKineticEnergy - postKineticEnergy;
          fEventAction->eDep += ed;
      }

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

          //fEventAction->numElasticSensitive += 1;
          //fEventAction->scatteredElastically = true;
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

       if (inelastic.compare(processName) == 0 && volumeName.compare(Absorber) == 0) {
           fEventAction->extInelastic = true;
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

