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
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include <sstream>
#include <unordered_map>
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
 fTrakLenCharged = fTrakLenNeutral = tZero = tOne = tHelp = eDep = enDep = tUp = tDown = coincidenceTime = 0.; 
 fNbStepsCharged = fNbStepsNeutral = numElastic = numElasticSensitive = detected = scatteredElastically = scatteredInelastically = nucleusRecoilEnergy = nScatterAngle = 0;
 nScint = nCher = numPhotons = largestTrackID = lastTracked = numReachedUp = numReachedDown = numPE = num = numSteps = totEnergy = 0;
 fTransmitFlag   = fReflectFlag    = scatteredNotSensitive = extInelastic = fail = numInelastic = nnPrime = 0;
 tZeroEx = tOneEx = 10000;
 detector = ""; 
 Coincedence = true;

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
  //if (tZeroEx != 10000){
      //std::cout << std::endl;
      //std::cout << "end of event, tZeroEx = " << tZeroEx << " ns" << std::endl;
      //std::cout << std::endl;
  //}
  const PrimaryGeneratorAction* prim = static_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const DetectorConstruction* det = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());


  if (enDep >= 400){
      detected = true;
      //if (tZeroEx != 10000){
          //std::cout << "evtID= " << eventID << ",detector = " << detector << ",tOneEx = " << tOneEx << ", tOnetrue = " << tOne << ", tZeroEx = " << tZeroEx << ", tZero true = " << tZero << std::endl;
      //}
      //std::cout << "event = " << eventID << ", detector = " << detector << ", TOFEx = " << tOneEx - tZeroEx << ", TOF truth = " << tOne - tZero << ", numPE = " << numReachedUp + numReachedDown << std::endl;
  }

  //std::stringstream ss;
  //ss << "numPE_" << prim->GetParticleGun()->GetParticleDefinition()->GetParticleName() << "_" << det->birks << "_" << det->absorptionLength[0]/CLHEP::cm << ".csv";
  //std::string filename = ss.str();
 
  //run->csvfi.open(filename,std::ios_base::app);
  //run->csvfi << eventID << "," << numReachedDown << "," << numReachedUp << "," << numReachedDown + numReachedUp << "," << std::endl;
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

  G4double earlyCutOff;
  G4double lateCutOff;

  //G4double tpcMinTOF = 58.5*CLHEP::ns;
  //G4double tpcMaxTOF = 62.*CLHEP::ns;

  std::vector<double> minTOFSND = {146 * CLHEP::ns, 146 * CLHEP::ns, 88 * CLHEP::ns, 76 * CLHEP::ns, 77 * CLHEP::ns, 77 * CLHEP::ns, 78 * CLHEP::ns, 79 * CLHEP::ns};
  std::vector<double> maxTOFSND = {155 * CLHEP::ns, 156 * CLHEP::ns, 93.5 * CLHEP::ns, 81 * CLHEP::ns, 82 * CLHEP::ns, 82 * CLHEP::ns, 83 * CLHEP::ns, 84 * CLHEP::ns};
    
  //for (double& num : minTOFSND) {
      //num += tpcMinTOF;
  //}

  //for (double& num : maxTOFSND){
      //num += tpcMaxTOF;
  //}

  std::unordered_map<std::string, int> detectorMap = {
        {"A0", 0}, {"A1", 1}, {"A2", 2}, {"A3", 3},
        {"A4", 4}, {"A5", 5}, {"A6", 6}, {"A7", 7}
    };

  if (detectorMap.find(detector) != detectorMap.end()) {
        int index = detectorMap[detector]; // Get the corresponding index

        // Assign early and late cut-off values based on the detector
        earlyCutOff = minTOFSND[index];
        lateCutOff = maxTOFSND[index];
  }
  //std::cout << "early cutoff = " << earlyCutOff << ", late cutoff = " << lateCutOff << std::endl;

  if (detected && Coincedence && ((tOne- tZero) > earlyCutOff) && ((tOne-tZero) < lateCutOff)){
      //std::cout << earlyCutOff << "," << lateCutOff << "," << detector << std::endl;
      run->csvfi.open("lightYield.csv",std::ios_base::app);
      run->csvfi << eventID << "," << numPhotons << "," << numReachedUp << "," << numReachedDown << "," << coincidenceTime << "," << numElasticSensitive << "," << scatteredInelastically << "," << scatteredNotSensitive << "," << nucleusRecoilEnergy << "," << tOneEx - tZeroEx << "," << detector << "," << std::endl;
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

