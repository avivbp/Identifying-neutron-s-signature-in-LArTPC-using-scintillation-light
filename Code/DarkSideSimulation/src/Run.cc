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
/// \file electromagnetic/TestEm5/src/Run.cc
/// \brief Implementation of the Run class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <sstream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), 
  fParticle(0), fEkin(0.)
{

  //csvfi.open("absTime.csv");
  //csvfi << "eventID" << "," << "tAbs" << "," << std::endl;
  //csvfi.close();

  //const PrimaryGeneratorAction* primaryGenAction
        //= static_cast<const PrimaryGeneratorAction*>(
        //    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 
  //std::stringstream ss;
  //ss << "numPE_" << primaryGenAction->GetParticleGun()->GetParticleDefinition()->GetParticleName() << "_" << det->birks << "_" << det->absorptionLength[0]/CLHEP::cm << ".csv";
  //std::string filename = ss.str();

  //csvfi.open(filename);
  //csvfi << "eventID" << "," << "numReachedDown" << "," << "numReachedUp" << "," << "numPE" << "," << std::endl;
  //csvfi.close();

  csvfi.open("photopeak.csv");
  csvfi << "eventID" << "," << "numReachedDown" << "," << "numReachedUp" << "," << "numPE" << "," << std::endl;
  csvfi.close();

  csvfi.open("energy_diff.csv");
  csvfi << "eventID" << "," << "energy" << "," << std::endl;
  csvfi.close();

  csvfi.open("phot_eKin.csv");
  csvfi << "electron_eKin" << "," << std::endl;
  csvfi.close();

  csvfi.open("elasticPhot.csv");
  csvfi << "nucleusKineticEnergykeV" << "," << "numPhotonsCreated" << "," << "nucleusEnergykeVPost" << "," << std::endl;
  csvfi.close();

  csvfi.open("photLoc.csv");
  csvfi << "eventID" << "," << "xPos" << "," << "yPos" << "," << "zPos" << "," << "tAbsorbed" << "," << std::endl;
  csvfi.close();

  csvfi.open("edep.csv");
  csvfi << "eventID" << "," << "eDep" << "," << std::endl;
  csvfi.close();

  csvfi.open("yPos.csv");
  csvfi << "eventID" << "," << "xPos" << "," << "yPos" << "," << "zPos" << "," << std::endl;
  csvfi.close();

  csvfi.open("singleElasticAngles.csv");
  csvfi << "eventID" << "," << "angle" << "," << std::endl;
  csvfi.close();

  csvfi.open("PMTS.csv");
  csvfi << "numPhotons" << "," << "numUp" << "," << "numDown" << "," << "eDep" << "," << std::endl;
  csvfi.close();

  csvfi.open("lightYield.csv");
  csvfi << "eventID" << "," << "numPhotons" << "," << "numUp" << "," << "numDown" << "," << "coincidenceTime" << "," << "numElasticSensitive" << "," << "scatteredInelastically" << "," << "scatteredNotSensitive" << "," << "nucleusRecoilEnergy" << "," << "timeOfFlight" << "," << std::endl;
  csvfi.close();

  csvfi.open("scintDist.csv");
  csvfi << "eventID" << "," << "numPhotons" << "," << "nScint" << "," << "sumStepLength" << "," << std::endl;
  csvfi.close();

  numInteractions = numElastic = numInelastic = numNNPrime = numNGamma = numNP = 0;

  csvfi.open("experimentalEnergies.csv");
  csvfi << "eventID" << "," << "experimentalRecoilEnergy" << "," << "timeOfFlight" << "," << "success" << "," << std::endl;
  csvfi.close();

  csfile.open("gammaEnergyDeposit.csv");
  csfile << "EventID" << "," << "totEnergyDeposit" << "secondaryGammasDeposit" << "," << std::endl;
  csfile.close();

  csvf.open("runStats.csv");
 // csvf << "neutronNum" << "," << "protonNum" << "," << "gammaNum" << "," << "Z" << "," << "A" << "," << std::endl;
 // csvf << "reaction" << "," << std::endl;
  csvf << "neutronNum" << "," << "protonNum" << "," << "sumStepLengths" << "," << "numInteractions" << "," << "numElastic" << "," << "primary" << "," << std::endl;
  csvf.close();
  
  csvfi.open("gammaEDist.csv");
  csvfi << "eventID" << "," << "Energy" << "," << "time" << "," << std::endl;
  csvfi.close();

  csvfi.open("neutronEscapeEnergy.csv");
  csvfi << "eventID" << "," << "neutronEscapeEnergy" << "," << "nPrime" << "," << std::endl;
  csvfi.close();

  csvfil.open("interactionStats.csv");
  csvfil << "distanceDiff" << "," << "timeDiff" << "," << std::endl;
  csvfil.close();

  csvfil.open("nPrimeEnergies.csv");
  csvfil << "eventID" << "," << "energy" << "," << std::endl;
  csvfil.close();

  csvfi.open("elasticEnergies.csv");
  csvfi << "eventID" << "," << "energy" << "," << "neutronAngle" << "," << "timeOfFlight" << "," << "external" << "," << "numElastic" << "," << std::endl;
  csvfi.close();

  csvfi.open("numElasticDetected.csv");
  csvfi << "eventID" << "," << "detected" << "," << "numElastic" << "," << std::endl;
  csvfi.close();

  csvfi.open("inelasticEnergies.csv");
  csvfi << "eventID" << "," << "energy" << "," << "neutronAngle" << "," << "timeOfFlight" << "," << std::endl;
  csvfi.close();

  csvfi.open("numInelasticDetected.csv");
  csvfi << "eventID" << "," << "detected" << "," << std::endl;
  csvfi.close();

  csvfi.open("photonEDIst.csv");
  csvfi << "eventID" << "," << "depositedEnergy" << "," << std::endl;
  csvfi.close();

  csvfi.open("eventTypes.csv");
  csvfi << "eventID" << "," << "type" << "," << std::endl;
  csvfi.close();

  fEnergyDeposit  = fEnergyDeposit2  = 0.;
  fTrakLenCharged = fTrakLenCharged2 = 0.;
  fTrakLenNeutral = fTrakLenNeutral2 = 0.;
  fNbStepsCharged = fNbStepsCharged2 = 0.;
  fNbStepsNeutral = fNbStepsNeutral2 = 0.;
  fMscProjecTheta = fMscProjecTheta2 = 0.;
  fMscThetaCentral = 0.;

  fNbGamma = fNbElect = fNbPosit = numSuc = numBackground = numPassed = dEdx = 0;
  single_elastic = multiple_elastic = elastic_NNPrime = multiple_NNPrime = no_interaction = numNNPrime = 0;
  fTransmit[0] = fTransmit[1] = fReflect[0] = fReflect[1] = 0;
  
  fMscEntryCentral = 0;
  
  fEnergyLeak[0] = fEnergyLeak[1] = fEnergyLeak2[0] = fEnergyLeak2[1] = 0.;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
  
  //compute theta0
  fMscThetaCentral = 3*ComputeMscHighland();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  fMscThetaCentral = localRun->fMscThetaCentral;

  // accumulate sums
  //
  fEnergyDeposit   += localRun->fEnergyDeposit;  
  fEnergyDeposit2  += localRun->fEnergyDeposit2;  
  fTrakLenCharged  += localRun->fTrakLenCharged;
  fTrakLenCharged2 += localRun->fTrakLenCharged2;   
  fTrakLenNeutral  += localRun->fTrakLenNeutral;  
  fTrakLenNeutral2 += localRun->fTrakLenNeutral2;
  fNbStepsCharged  += localRun->fNbStepsCharged;
  fNbStepsCharged2 += localRun->fNbStepsCharged2;
  fNbStepsNeutral  += localRun->fNbStepsNeutral;
  fNbStepsNeutral2 += localRun->fNbStepsNeutral2;
  fMscProjecTheta  += localRun->fMscProjecTheta;
  fMscProjecTheta2 += localRun->fMscProjecTheta2;

    
  fNbGamma += localRun->fNbGamma;
  fNbElect += localRun->fNbElect;      
  fNbPosit += localRun->fNbPosit;
  
  fTransmit[0] += localRun->fTransmit[0];  
  fTransmit[1] += localRun->fTransmit[1];
  fReflect[0]  += localRun->fReflect[0];
  fReflect[1]  += localRun->fReflect[1];
  
  fMscEntryCentral += localRun->fMscEntryCentral;
  
  fEnergyLeak[0]  += localRun->fEnergyLeak[0];
  fEnergyLeak[1]  += localRun->fEnergyLeak[1];
  fEnergyLeak2[0] += localRun->fEnergyLeak2[0];
  fEnergyLeak2[1] += localRun->fEnergyLeak2[1];
                  
  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  csvf.close();
  csvfi.close();
  csvfil.close();
 // std::cout << "total number of interactions: " << numInteractions << std::endl;
  //std::cout << "number of elastic interactions: " << numElastic << std::endl;
  //std::cout << "number of inelastic interactions: " << numInelastic << std::endl;
  std::cout << "number of optical photons passed = " << numPassed << std::endl;
  //std::cout << "number of single elastic events with 45-55 keV Ar recoil nucleus is : " << numSuc << std::endl;
  //std::cout << "number of background events that entered the liquid scintillator is : " << numBackground << std::endl;

  //std::cout << "number of single elastic events = " << single_elastic << std::endl;
  //std::cout << "number of multiple elastic events = " << multiple_elastic << std::endl;
  //std::cout << "number of multiple elastic + (n,n') events = " << elastic_NNPrime << std::endl;
  //std::cout << "number of multiple (n,n') events = " << multiple_NNPrime << std::endl;
  //std::cout << "number of events with no interaction = " << no_interaction << std::endl;
  // compute mean and rms
  //
  G4int TotNbofEvents = numberOfEvent;
  if (TotNbofEvents == 0) return;
  
  G4double EnergyBalance = fEnergyDeposit + fEnergyLeak[0] + fEnergyLeak[1];
  EnergyBalance /= TotNbofEvents;

  fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit*fEnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep/TotNbofEvents);
  else            rmsEdep = 0.;

  fTrakLenCharged /= TotNbofEvents; fTrakLenCharged2 /= TotNbofEvents;
  G4double rmsTLCh = fTrakLenCharged2 - fTrakLenCharged*fTrakLenCharged;
  if (rmsTLCh>0.) rmsTLCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsTLCh = 0.;

  fTrakLenNeutral /= TotNbofEvents; fTrakLenNeutral2 /= TotNbofEvents;
  G4double rmsTLNe = fTrakLenNeutral2 - fTrakLenNeutral*fTrakLenNeutral;
  if (rmsTLNe>0.) rmsTLNe = std::sqrt(rmsTLNe/TotNbofEvents);
  else            rmsTLNe = 0.;

  fNbStepsCharged /= TotNbofEvents; fNbStepsCharged2 /= TotNbofEvents;
  G4double rmsStCh = fNbStepsCharged2 - fNbStepsCharged*fNbStepsCharged;
  if (rmsStCh>0.) rmsStCh = std::sqrt(rmsStCh/TotNbofEvents);
  else            rmsStCh = 0.;

  fNbStepsNeutral /= TotNbofEvents; fNbStepsNeutral2 /= TotNbofEvents;
  G4double rmsStNe = fNbStepsNeutral2 - fNbStepsNeutral*fNbStepsNeutral;
  if (rmsStNe>0.) rmsStNe = std::sqrt(rmsStNe/TotNbofEvents);
  else            rmsStNe = 0.;

  G4double Gamma = (G4double)fNbGamma/TotNbofEvents;
  G4double Elect = (G4double)fNbElect/TotNbofEvents;
  G4double Posit = (G4double)fNbPosit/TotNbofEvents;

  G4double transmit[2];
  transmit[0] = 100.*fTransmit[0]/TotNbofEvents;
  transmit[1] = 100.*fTransmit[1]/TotNbofEvents;

  G4double reflect[2];
  reflect[0] = 100.*fReflect[0]/TotNbofEvents;
  reflect[1] = 100.*fReflect[1]/TotNbofEvents;

  G4double rmsMsc = 0., tailMsc = 0.;
  if (fMscEntryCentral > 0) {
    fMscProjecTheta /= fMscEntryCentral; fMscProjecTheta2 /= fMscEntryCentral;
    rmsMsc = fMscProjecTheta2 - fMscProjecTheta*fMscProjecTheta;
    if (rmsMsc > 0.) { rmsMsc = std::sqrt(rmsMsc); }
    if(fTransmit[1] > 0.0) {
      tailMsc = 100.- (100.*fMscEntryCentral)/(2*fTransmit[1]);
    }    
  }
  
  fEnergyLeak[0] /= TotNbofEvents; fEnergyLeak2[0] /= TotNbofEvents;
  G4double rmsEl0 = fEnergyLeak2[0] - fEnergyLeak[0]*fEnergyLeak[0];
  if (rmsEl0>0.) rmsEl0 = std::sqrt(rmsEl0/TotNbofEvents);
  else           rmsEl0 = 0.;
  
  fEnergyLeak[1] /= TotNbofEvents; fEnergyLeak2[1] /= TotNbofEvents;
  G4double rmsEl1 = fEnergyLeak2[1] - fEnergyLeak[1]*fEnergyLeak[1];
  if (rmsEl1>0.) rmsEl1 = std::sqrt(rmsEl1/TotNbofEvents);
  else           rmsEl1 = 0.;    
  
      
  //Stopping Power from input Table.
  //
  const G4Material* material = fDetector->GetAbsorberMaterial();
  G4double length      = fDetector->GetAbsorberThickness();
  G4double density     = material->GetDensity();   
  G4String partName    = fParticle->GetParticleName();

  G4EmCalculator emCalculator;
  G4double dEdxTable = 0., dEdxFull = 0.;
  if (fParticle->GetPDGCharge()!= 0.) { 
    dEdxTable = emCalculator.GetDEDX(fEkin,fParticle,material);
    dEdxFull  = emCalculator.ComputeTotalDEDX(fEkin,fParticle,material);    
  }
  G4double stopTable = dEdxTable/density;
  G4double stopFull  = dEdxFull /density; 
   
  //Stopping Power from simulation.
  //    
  G4double meandEdx  = fEnergyDeposit/length;
  G4double stopPower = meandEdx/density;  

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(3);
  
  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  G4cout.precision(4);
  
  G4cout << "\n Total energy deposit in absorber per event = "
         << G4BestUnit(fEnergyDeposit,"Energy") << " +- "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
         
  G4cout << "\n -----> Mean dE/dx = " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "\n From formulas :" << G4endl; 
  G4cout << "   restricted dEdx = " << dEdxTable/(MeV/cm) << " MeV/cm"
         << "\t(" << stopTable/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "   full dEdx       = " << dEdxFull/(MeV/cm) << " MeV/cm"
         << "\t(" << stopFull/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;
         
  G4cout << "\n Leakage :  primary = "
         << G4BestUnit(fEnergyLeak[0],"Energy") << " +- "
         << G4BestUnit(rmsEl0,       "Energy")
         << "   secondaries = "
         << G4BestUnit(fEnergyLeak[1],"Energy") << " +- "
         << G4BestUnit(rmsEl1,       "Energy")          
         << G4endl;
         
  G4cout << " Energy balance :  edep + eleak = "
         << G4BestUnit(EnergyBalance,"Energy")
         << G4endl;         
                           
  G4cout << "\n Total track length (charged) in absorber per event = "
         << G4BestUnit(fTrakLenCharged,"Length") << " +- "
         << G4BestUnit(rmsTLCh,       "Length") << G4endl;

  G4cout << " Total track length (neutral) in absorber per event = "
         << G4BestUnit(fTrakLenNeutral,"Length") << " +- "
         << G4BestUnit(rmsTLNe,       "Length") << G4endl;

  G4cout << "\n Number of steps (charged) in absorber per event = "
         << fNbStepsCharged << " +- " << rmsStCh << G4endl;

  G4cout << " Number of steps (neutral) in absorber per event = "
         << fNbStepsNeutral << " +- " << rmsStNe << G4endl;

  G4cout << "\n Number of secondaries per event : Gammas = " << Gamma
         << ";   electrons = " << Elect
           << ";   positrons = " << Posit << G4endl;

  G4cout << "\n Number of events with the primary particle transmitted = "
         << transmit[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle transmitted "
         << "(same charge as primary) = " << transmit[0] << " %" << G4endl;

  G4cout << "\n Number of events with the primary particle reflected = "
         << reflect[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle reflected "
         << "(same charge as primary) = " << reflect[0] << " %" << G4endl;

  // compute width of the Gaussian central part of the MultipleScattering
  //
  G4cout << "\n MultipleScattering:" 
         << "\n  rms proj angle of transmit primary particle = "
         << rmsMsc/mrad << " mrad (central part only)" << G4endl;

  G4cout << "  computed theta0 (Highland formula)          = "
         << ComputeMscHighland()/mrad << " mrad" << G4endl;
           
  G4cout << "  central part defined as +- "
         << fMscThetaCentral/mrad << " mrad; " 
         << "  Tail ratio = " << tailMsc << " %" << G4endl;
         
  // normalize histograms
  //
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  //G4int ih = 1;
  //G4double binWidth = analysisManager->GetH1Width(ih);
  //G4double fac = 1./(TotNbofEvents*binWidth);
 // analysisManager->ScaleH1(ih,fac);

  //ih = 10;
  //binWidth = analysisManager->GetH1Width(ih);
  //fac = 1./(TotNbofEvents*binWidth);
  //analysisManager->ScaleH1(ih,fac);

  //ih = 12;
  //analysisManager->ScaleH1(ih,1./TotNbofEvents);

  // reset default precision
  G4cout.precision(prec);
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Run::ComputeMscHighland()
{
  //compute the width of the Gaussian central part of the MultipleScattering
  //projected angular distribution.
  //Eur. Phys. Jour. C15 (2000) page 166, formule 23.9

  G4double t = (fDetector->GetAbsorberThickness())
              /(fDetector->GetAbsorberMaterial()->GetRadlen());
  if (t < DBL_MIN) return 0.;

  G4double T = fEkin;
  G4double M = fParticle->GetPDGMass();
  G4double z = std::abs(fParticle->GetPDGCharge()/eplus);

  G4double bpc = T*(T+2*M)/(T+M);
  G4double teta0 = 13.6*MeV*z*std::sqrt(t)*(1.+0.038*std::log(t))/bpc;
  return teta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......