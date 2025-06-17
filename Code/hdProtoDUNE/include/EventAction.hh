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
/// \file electromagnetic/TestEm5/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <list>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    G4int numInteractions;
    G4bool wroteToFile;
    G4double sumStepLengths;
    G4int totNeutronNum;
    G4int totProtonNum;
    G4double gammaEnergyDeposit;
    G4double secondarygammasDeposit;
    G4bool passed;
    G4bool performedNPrime;
    G4int numElastic;
    G4bool detected;
    G4bool scatteredElastically;
    G4int numElasticSensitive;
    G4bool scatteredInelastically;
    G4double nucleusRecoilEnergy;
    G4double nScatterAngle;
    G4double tZero;
    G4double tOne;
    G4bool scatteredNotSensitive;
    G4bool extInelastic;
    G4bool Coincedence;
    G4bool fail;
    G4bool aborted;
    G4bool nCapture;
    G4bool secondaryNeutron;
    CLHEP::Hep3Vector pos0;
    CLHEP::Hep3Vector pos1;
    G4double experimentalRecoilEnergy;
    G4double tZeroEx;
    G4double tOneEx;
    G4double tHelp;
    G4int nScint;
    G4int nCher;
    G4int numPhotons;
    G4int numPE;
    G4int largestTrackID;
    G4int lastTracked;
    G4int numReachedUp;
    G4int numReachedDown;
    G4int num;
    G4int numInelastic;
    G4int nnPrime;
    G4int numSteps;
    G4int numSurface;
    G4int numBases;
    G4double totEnergy;
    G4double tUp;
    G4double tDown;
    G4double eDep;
    G4double enDep;
    G4double sumStepLength;
    G4double coincidenceTime;
    G4String detector;
    std::ofstream csvfile;
    std::list<int> parentIDS;
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
 
    void numInteractionsPP(){numInteractions += 1;};
    void setNPrime(G4bool b){performedNPrime = b;};
    void resetNumInteractions(){numInteractions = 0;};
    void setPassed(G4bool b){passed = b;};
    void addStepLength(G4double stepL){sumStepLengths+=stepL;};
    void addNeutronNum(G4int num){totNeutronNum += num;};
    void addProtonNum(G4int num){totProtonNum += num;};
    void resetStepLengths(){sumStepLengths = 0.;};
    void resetNeutronNum(){totNeutronNum = 0;};
    void resetProtonNum(){totProtonNum = 0;};
    void setWroteToFile(G4bool boo ){wroteToFile = boo;};
    void addEDep(G4double gammaEDep){gammaEnergyDeposit+=gammaEDep;};
    void addDeposit(G4double gammaEDep){secondarygammasDeposit+=gammaEDep;};
    void elasticPP(){numElastic += 1;};

    void AddEnergy      (G4double edep)   {fEnergyDeposit  += edep;};
    void AddTrakLenCharg(G4double length) {fTrakLenCharged += length;};
    void AddTrakLenNeutr(G4double length) {fTrakLenNeutral += length;};
    
    void CountStepsCharg ()               {fNbStepsCharged++ ;};
    void CountStepsNeutr ()               {fNbStepsNeutral++ ;};
    
    void SetTransmitFlag (G4int flag) 
                           {if (flag > fTransmitFlag) fTransmitFlag = flag;};
    void SetReflectFlag  (G4int flag) 
                           {if (flag > fReflectFlag)   fReflectFlag = flag;};
                                             
        
  private:
    G4double fEnergyDeposit;
    G4double fTrakLenCharged, fTrakLenNeutral;
    G4int    fNbStepsCharged, fNbStepsNeutral;
    G4int    fTransmitFlag,   fReflectFlag;        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
