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
/// \file electromagnetic/TestEm5/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
//G4bool wroteDistance = false;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
 fEnergyDeposit(0.),
 fTrakLenCharged(0.), fTrakLenNeutral(0.),
 fNbStepsCharged(0), fNbStepsNeutral(0),
 fTransmitFlag(0), fReflectFlag(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt )
{
// std::cout << "in begin of event action" << std::endl;
 // initialisation per event
 setWroteToFile(false);
 resetStepLengths();
 resetNumInteractions();
 resetNeutronNum();
 resetProtonNum();
 setNPrime(false);
 gammaEnergyDeposit = secondarygammasDeposit = 0.;
// std::cout << passed << std::endl;
// setPassed(false);
// std::cout << passed << std::endl;
 fEnergyDeposit  = experimentalRecoilEnergy = sumStepLength = 0.;
 pos0 = pos1 = G4ThreeVector(0.,0.,0.);
 fTrakLenCharged = fTrakLenNeutral = tZero = tOne = tZeroEx = tHelp = eDep = enDep = tUp = tDown = coincidenceTime = 0.; 
 fNbStepsCharged = fNbStepsNeutral = numElastic = numElasticSensitive = detected = scatteredElastically = scatteredInelastically = nucleusRecoilEnergy = nScatterAngle = 0;
 nScint = nCher = numPhotons = largestTrackID = lastTracked = numReachedUp = numReachedDown = numPE = num = numSteps = totEnergy = 0;
 fTransmitFlag   = fReflectFlag    = scatteredNotSensitive = extInelastic = Coincedence = fail = numInelastic = nnPrime = 0;

 csvfile.open("scintCheck.csv");
 csvfile << "energy" << "," << "time" << "," << "theta" << "," << "phi" << "," << "xpos" << "," << "ypos" << "," << "zpos" << "," << std::endl;
 csvfile.close();

 csvfile.open("primaryPos.csv");
 csvfile << "eventID" << "," << "xpos" << "," << "ypos" << "," << "zpos" << "," << std::endl;
 csvfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int eventID = evt->GetEventID();

  //std::cout << "num photons reached top PMT array = " << numReachedUp << std::endl;
  //std::cout << "num photons reached bottom PMT = " << numReachedDown << std::endl;
  //std::cout << "tot num = " << numReachedDown + numReachedUp << std::endl;

  if (eDep < 550 && eDep > 450){
      run->csvfi.open("photopeak.csv",std::ios_base::app);
      run->csvfi << "eventID" << "," << "numReachedDown" << "," << "numReachedUp" << "," << "numPE" << "," << std::endl;
      run->csvfi.close();
  }

  const DetectorConstruction* det = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  const PrimaryGeneratorAction* prim = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  std::stringstream ss;
  ss << "numPE_" << prim->GetParticleGun()->GetParticleDefinition()->GetParticleName() << "_" << det->birks << "_" << det->absorptionLength[0]/CLHEP::cm << ".csv";
  std::string filename = ss.str();

  if (eventID == 0){
      run->csvfi.open(filename);
      run->csvfi << "eventID" << "," << "numReachedDown" << "," << "numReachedUp" << "," << "numPE" << "," << std::endl;
      run->csvfi.close();
  }

  run->csvfi.open(filename,std::ios_base::app);
  run->csvfi << eventID << "," << numReachedDown << "," << numReachedUp << "," << numReachedDown + numReachedUp << "," << std::endl;
  run->csvfi.close();

  if (!fail && numElasticSensitive == 1){
      //std::cout << "single elastic 45-55 keV" << std::endl;
      if (detected){
          run->numSuc += 1;
      }
      //std::cout << "nice" << std::endl;
  }

  if (fail && detected && !scatteredInelastically && numElasticSensitive == 1){
      run->numBackground += 1;
  }

  //run->csvfi.open("eventTypes.csv",std::ios_base::app);

  //if (numElastic == 1 && numInelastic == 0){
      //run->single_elastic +=1;
      //run->csvfi << eventID << "," << "single_elastic" << "," << std::endl;     
      //run->csvfi.open("singleElasticAngles.csv",std::ios_base::app);
      //run->csvfi << evt->GetEventID() << "," << nScatterAngle << "," << std::endl;
      //run->csvfi.close();
  //}
  //else if (numElastic > 1 && numInelastic == 0){
      //run->multiple_elastic += 1;
      //run->csvfi << eventID << "," << "multiple_elastic" << "," << std::endl;
  //}
  //else if (numElastic >= 1 && numInelastic >= 1){
      //run->elastic_NNPrime += 1;
      //run->csvfi << eventID << "," << "elastics_plus_nnprime" << "," << std::endl;
  //}
  //else if (numInelastic > 1){
      //run->multiple_NNPrime += 1;
      //run->csvfi << eventID << "," << "multiple_nnprime" << "," << std::endl;
  //}

  //else if (numInelastic == 1 && numElastic == 0){
      //run-> numNNPrime += 1;
      //run->csvfi << eventID << "," << "nnPrime" << "," << std::endl;
  //}

  //else {
      //run->no_interaction += 1;
      //run->csvfi << eventID << "," << "no interaction" << "," << std::endl;
  //}
  //run->csvfi.close();

 // std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout<<std::endl;
  //std::cout << "number of optical photons created by lower energy electrons from phot interaction = " << num << std::endl;
  //std::cout << "number of photons that reached top PMT = " << numReachedUp << std::endl;
  //std::cout << "number of photons that reached bottom PMT = " << numReachedDown << std::endl;
  //std::cout << "energy deposit by electrons = " << eDep << " MeV" << std::endl;
  //std::cout << "number of scintillation photons created = " << numPhotons << std::endl;
  //std::cout << "number of p.e = " << numPE << std::endl;
  //std::cout << "number of cherenkov photons = " << nCher << std::endl;
  //std::cout << "sumStepLength = " << sumStepLength << std::endl;
  //std::cout << std::endl;

 
  if (detected && Coincedence && ((tOne - tZero) > 42.*CLHEP::ns) && ((tOne-tZero) < 49.*CLHEP::ns)){
      run->csvfi.open("lightYield.csv",std::ios_base::app);
      run->csvfi << eventID << "," << numPhotons << "," << numReachedUp << "," << numReachedDown << "," << coincidenceTime << "," << numElasticSensitive << "," << scatteredInelastically << "," << scatteredNotSensitive << "," << nucleusRecoilEnergy << "," << tOne - tZero << "," << std::endl;
      run->csvfi.close();
  }

  if (eventID%10000 == 0){
      std::cout << "event number " << evt->GetEventID() << std::endl;
  }


 run->AddEnergy(fEnergyDeposit);
 run->AddTrakLenCharg(fTrakLenCharged);
 run->AddTrakLenNeutr(fTrakLenNeutral);

 run->CountStepsCharg(fNbStepsCharged);
 run->CountStepsNeutr(fNbStepsNeutral);

 run->CountTransmit (fTransmitFlag);
 run->CountReflect  (fReflectFlag);
 
// if (fEnergyDeposit > 0.){
 //   G4AnalysisManager::Instance()->FillH1(1,fEnergyDeposit);
  //  G4AnalysisManager::Instance()->FillNtupleDColumn(1,fEnergyDeposit);
   // G4AnalysisManager::Instance()->AddNtupleRow();
 //   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

