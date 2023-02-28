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

#include "Run.hh"
#include "HistoManager.hh"

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
 fTrakLenCharged = fTrakLenNeutral = tZero = tOne = tZeroEx = tHelp = eDep = tUp = tDown = coincidenceTime = 0.; 
 fNbStepsCharged = fNbStepsNeutral = numElastic = numElasticSensitive = detected = scatteredElastically = scatteredInelastically = nucleusRecoilEnergy = nScatterAngle = 0;
 nScint = nCher = numPhotons = largestTrackID = lastTracked = numReachedUp = numReachedDown = 0;
 fTransmitFlag   = fReflectFlag    = scatteredNotSensitive = extInelastic = Coincedence = 0;    

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
 // std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
  //std::cout<<std::endl;

  //std::cout << "number of photons that reached top PMT = " << numReachedUp << std::endl;
  //std::cout << "number of photons that reached bottom PMT = " << numReachedDown << std::endl;
  //std::cout << "energy deposit by electrons = " << eDep << " MeV" << std::endl;
  //std::cout << "number of scintillation photons by hand = " << numPhotons << std::endl;
  //std::cout << "number of cherenkov photons = " << nCher << std::endl;
  //std::cout << "sumStepLength = " << sumStepLength << std::endl;
  //std::cout << std::endl;

  if (!wroteToFile){
      //csvfile.open("neutronEscapeEnergy.csv",std::ios_base::app);
      //csvfile << evt->GetEventID() << "," << 0 << "," << 0 << "," << std::endl;
      if (performedNPrime){
          //csvfile << evt->GetEventID() << "," << 0 << "," << 1 << "," << std::endl;
      }
      //csvfile.close();
  }
 
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int eventID = evt->GetEventID();

  //run->csvfi.open("scintDist.csv",std::ios_base::app);
  //run->csvfi << eventID << "," << numPhotons << "," << nScint << "," << sumStepLength << "," << std::endl;
  //run->csvfi.close();

  run->csvfi.open("PMTS.csv",std::ios_base::app);
  run->csvfi << numPhotons << "," << numReachedUp << "," << numReachedDown << "," << eDep << "," << std::endl;
  run->csvfi.close();

  if (detected && Coincedence && ((tOne - tZero) > 42.*CLHEP::ns) && ((tOne-tZero) < 49.*CLHEP::ns)){
      run->csvfi.open("lightYield.csv",std::ios_base::app);
      run->csvfi << eventID << "," << numPhotons << "," << numReachedUp << "," << numReachedDown << "," << coincidenceTime << "," << numElasticSensitive << "," << scatteredInelastically << "," << scatteredNotSensitive << "," << std::endl;
      run->csvfi.close();
  }

  // only if there is energy deposit in the sensitive volume will there be scintillation light
  //if (detected && (numElasticSensitive >= 1 || scatteredInelastically)){
      //run->csvfi.open("experimentalEnergies.csv",std::ios_base::app);
      //run->csvfi << eventID << "," << experimentalRecoilEnergy << "," << tOne - tZeroEx << "," << std::endl;
      //run->csvfi.close();
  //}

  if (eventID%100000 == 0){
      std::cout << "event number " << evt->GetEventID() << std::endl;
  }

  //if ( (nucleusRecoilEnergy < 5 || nucleusRecoilEnergy > 20) && detected && numElasticSensitive == 1 && !scatteredInelastically && !extInelastic){
      //std::cout << "sus event : " << eventID << std::endl;
  //}

  //if (detected && !scatteredInelastically) {

      //if ((tOne - tZero) < 40){
          //std::cout << " too fast in eventID = " << eventID << std::endl;
      //}

      //run->csvfi.open("elasticEnergies.csv",std::ios_base::app);
      //run->csvfi << evt->GetEventID() << "," << nucleusRecoilEnergy << "," << nScatterAngle << "," << tOne - tZero  << "," << scatteredNotSensitive << "," << numElasticSensitive << "," << std::endl;
      //run->csvfi.close();
  //}

  //if (detected && scatteredInelastically){
    
      //run->csvfi.open("inelasticEnergies.csv",std::ios_base::app);
      //run->csvfi << evt->GetEventID() << "," << nucleusRecoilEnergy << "," << nScatterAngle << "," << tOne - tZero << "," << std::endl;
      //run->csvfi.close();     
  //}

  //csvfile.open("runStats.csv",std::ios_base::app);
  //csvfile << totNeutronNum << "," << totProtonNum << "," << sumStepLengths << "," << numInteractions << "," << numElastic << "," << !performedNPrime << "," << std::endl;
  //csvfile.close();

  //csvfile.open("gammaEnergyDeposit.csv",std::ios_base::app);
  //csvfile << evt->GetEventID() << "," << gammaEnergyDeposit << "," << secondarygammasDeposit << "," << std::endl;
  //csvfile.close();

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

