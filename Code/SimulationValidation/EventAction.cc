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

void EventAction::BeginOfEventAction(const G4Event* )
{
// std::cout << "in begin of event action" << std::endl;
 // initialisation per event
 setWroteToFile(false);
 resetStepLengths();
 resetNumInteractions();
// std::cout << passed << std::endl;
// setPassed(false);
// std::cout << passed << std::endl;
 fEnergyDeposit  = 0.;
 fTrakLenCharged = fTrakLenNeutral = 0.; 
 fNbStepsCharged = fNbStepsNeutral = 0;
 fTransmitFlag   = fReflectFlag    = 0;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 // csvfile.open("results.csv",std::ios_base::app);
  //csvfile << evt->GetEventID() << "," << fEnergyDeposit << "," << std::endl;
 // csvfile.open("freePathDist.csv",std::ios_base::app);
  //csvfile << evt->GetEventID() << "," << 0 << "," << std::endl;
  //csvfile.close();
 // std::cout << std::endl;
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
// if (passed){
 //    run->addNumPassed();
// }             
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

